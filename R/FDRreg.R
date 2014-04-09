FDRreg = function(z, features, nulltype='empirical', nmids=150, densknots=10, MCMC=NULL, center=TRUE, scale=TRUE) {
# False discovery rate regression
# z = vector of z scores
# X = design matrix of covariates, assumed NOT to have an intercept just as in vanilla lm()
# nulltype = flag for what kind of null hypothesis to assume, theoretical or empirical
	N = length(z)
	
	X = cbind(1,scale(features, center=center, scale=scale))
	P = ncol(X)
	
	# Estimate the marginal density
	l1 = efron(z, nmids=nmids, df=densknots, nulltype=nulltype)
	mu0 = l1$mu0
	sig0 = l1$sig0
	p0 = l1$p0
	M0 = dnorm(z, mu0, sig0)
	grid.f0z = dnorm(l1$mids,mu0,sig0)
	MTot = l1$fz
	nullcases = which(p0*M0>MTot)
	
	# Iteratively refine estimate for beta
	betaguess = rep(0,P)
	travel=1
	W = (1-p0)
	while(travel > 1e-6) {
		PostProb = pmax(0, 1 - (1-W)*M0/MTot)
		PostProb[nullcases] = 0
		suppressWarnings(lm1 <- glm(PostProb ~ X-1, family=binomial))
		newbeta = coef(lm1)
		W = ilogit(X %*% newbeta)
		travel = sum(abs(betaguess-newbeta))
		betaguess = newbeta
	}
	PostProb = pmax(0, 1 - (1-W)*M0/MTot)
	PostProb[nullcases] = 0
	out2 = getFDR(PostProb)
	
	if(missing(MCMC)) {
	  betasave=NULL
	}
	else {
	  # Initialize MCMC
	  PriorPrec = diag(c(1e-9, rep(0.01, P-1)))
	  PriorMean = rep(0,P)
	
	  # Pass to C++
	  out1 = FDRregCPP(z, X, M0=M0, MTot = MTot, PriorPrec, PriorMean, nmc=MCMC$nmc, nburn=MCMC$burn, p0=p0, betaguess = betaguess)
	  PostProb = out1$postprob
	  out2 = getFDR(PostProb)
	  betasave = out1$betasave
	}

	list(z=z, localfdr=out2$localfdr, FDR=out2$FDR, X=X, grid=l1$mids, breaks=l1$breaks,
         grid.fz=l1$zdens, grid.f0z=grid.f0z, grid.zcounts=l1$zcounts, dnull = M0, dmix=MTot,
         mu0=mu0, sig0=sig0, p0=p0, priorprob = W, postprob = PostProb,
         model=lm1, betasave = betasave
    )
}



BayesFDRreg = function(z, features, nulltype='empirical', nmc=10000, nburn=500, nmids=150,
	densknots=10, ncomps=NULL, priorpars=NULL, fullbayes = FALSE, center=TRUE, scale=TRUE, verbose=NULL) {
# More Bayesian version of false discovery rate regression
# z = vector of z scores
# features = design matrix of covariates, assumed NOT to have an intercept just as in vanilla lm()
# nulltype = flag for what kind of null hypothesis to assume, theoretical or empirical
	N = length(z)
	X = cbind(1,scale(features, center=center, scale=scale))
	P = ncol(X)
	
	if(missing(verbose)) verbose = nmc+nburn+1
	
	# Estimate the marginal density
	l1 = efron(z, nmids=nmids, df=densknots, nulltype=nulltype)
	mu0 = l1$mu0
	sig0 = l1$sig0
	p0 = l1$p0
	M0 = dnorm(z, mu0, sig0)
	grid.f0z = dnorm(l1$mids,mu0,sig0)	
	MTot = l1$fz
		
	# Initialize MCMC
	if(missing(priorpars)) {
		PriorPrec = diag(c(1e-9, rep(1/9, P-1)))
		PriorMean = rep(0,P)
	} else{
		PriorPrec = priorpars$PriorPrec
		PriorMean = priorpars$PriorMean
	}
		
	if(!fullbayes) {
		
		# Use EM to fit the alternative density
		if(missing(ncomps)) {
			foundfit = FALSE
			ncomps = 1
			emfit = deconvolveEM(z, ncomps)
			while(!foundfit) {
				newfit = deconvolveEM(z, ncomps+1)
				if(newfit$AIC > emfit$AIC) {
					foundfit = TRUE
				} else {
					emfit = newfit
					ncomps = ncomps+1
				}
			}
		} else {
			emfit = deconvolveEM(z, ncomps)
		}
		M1 = dnormix(z, emfit$weights[1:ncomps]/sum(emfit$weights[1:ncomps]), emfit$means[1:ncomps], emfit$vars[1:ncomps])
	
		postprob_init = (1-p0)*M1/{(1-p0)*M1 + p0*M0}
		suppressWarnings(lm1 <- glm(postprob_init ~ X-1, family=binomial))
		
		# Pass to C++ and return result
		out1 = EmpiricalBayesFDRregCPP(z=z, X=X, M0=M0, M1=M1,PriorPrecision=PriorPrec, PriorMean=PriorMean,
			nmc=nmc, nburn=nburn, betaguess = coef(lm1))
		out2 = getFDR(out1$postprob)
		mylist = list(z=z, localfdr=out2$localfdr, FDR=out2$FDR, X=X, grid=l1$mids, breaks=l1$breaks,
        	grid.fz=l1$zdens, grid.f0z=grid.f0z, grid.zcounts=l1$zcounts,
        	dnull = M0, dmix=MTot, dalt = M1, mu0=mu0, sig0=sig0, p0=p0, ncomps=ncomps,
        	betasave = out1$betasave,
        	priorprob = out1$priorprob, postprob = out1$postprob
        )
    }

	else {
		
		if(missing(ncomps)) {
			foundfit = FALSE
			ncomps = 1
			emfit = deconvolveEM(z, ncomps)
			while(!foundfit) {
				newfit = deconvolveEM(z, ncomps+1)
				if(newfit$AIC > emfit$AIC) {
					foundfit = TRUE
				} else {
					emfit = newfit
					ncomps = ncomps+1
				}
			}
			M1 = dnormix(z, emfit$weights[1:ncomps]/sum(emfit$weights[1:ncomps]), emfit$means[1:ncomps], emfit$vars[1:ncomps])
		} else 	M1 = dnorm(z, 0, 4)
		
		PriorPrecXMean = PriorPrec %*% PriorMean
		Beta = rep(0,P)
		Beta[1] = flogit(sum(abs(z)>2)/N)
		BetaSave = matrix(0, nrow=nmc, ncol=P)
		MuSave = matrix(0, nrow=nmc, ncol=ncomps)
		VarSave = matrix(0, nrow=nmc, ncol=ncomps)
		WeightsSave = matrix(0, nrow=nmc, ncol=ncomps)
		PostProbSave = 0
		PriorProbSave = 0
		M1Save = 0
		
		# Alternative hypothesis
		comp_weights = rep(1/ncomps, ncomps)
		comp_means = quantile(z[abs(z)>2], probs=seq(0.025,0.975,length=ncomps))
		comp_variance = rep(1, ncomps)
		myvar = comp_variance + 1
		
		# Main MCMC
		for(t in 1:(nmc+nburn)) {
			
			if(t %% verbose == 0) cat(t, "\n")
			
			### Update indicators
			Psi = drop(X %*% Beta)
			W = ilogit(Psi)
			PostProb = W*M1/{(1-W)*M0 + W*M1}		
			Gamma = rbinom(N,1,PostProb)
			
			
			### Update mixture of normals model
			
			# Sufficient statistics
			signals = z[which(Gamma==1)]
			components = draw_mixture_component(signals, rep(1,N), weights=comp_weights, mu = comp_means, tau2 = comp_variance) + 1
			sumsignals = mosaic::maggregate(signals ~ factor(components, levels=1:ncomps), FUN='sum')
			nsig = mosaic::maggregate(signals ~ factor(components, levels=1:ncomps), FUN='length')
			tss = mosaic::maggregate((signals-comp_means[components])^2 ~ factor(components, levels=1:ncomps), FUN='sum')
			
			#print(c(comp_weights, comp_means, comp_variance))
			
			# Updates
			for(k in 1:ncomps) myvar[k] = 1/rtgamma_once({nsig[k]+0.5}/2, {tss[k]+0.5}/2, 0, 1)
			comp_variance = myvar - 1
			muvar = comp_variance/{nsig + comp_variance*0.01}
			muhat = sumsignals/{nsig + comp_variance*0.01}
			comp_means = rnorm(ncomps, muhat, sqrt(muvar))	
			comp_weights = rdirichlet_once(rep(1, ncomps) + nsig)
			M1 = dnormix(z, comp_weights, comp_means, myvar)
			
			### Update latent variables in logit likelihood
			Om = as.numeric(BayesLogit::rpg(N,rep(1,N),Psi))
	
			### Update regression parameters
			Kap = PostProb - 1/2
			PrecMat = t(X) %*% {Om * X} + PriorPrec
			Beta.V = solve(PrecMat)
			Beta.mu = Beta.V %*% {t(X) %*% Kap + PriorPrecXMean}
			Beta = t(mvtnorm::rmvnorm(1,mean=Beta.mu,sigma=Beta.V))	
			if(t > nburn) {
				BetaSave[t-nburn,] = Beta
				MuSave[t-nburn,] = comp_means
				VarSave[t-nburn,] = comp_variance
				WeightsSave[t-nburn,] = comp_weights
				PostProbSave = PostProbSave + (1/nmc)*PostProb
				PriorProbSave = PriorProbSave + (1/nmc)*W
				M1Save = M1Save + (1/nmc)*M1
			}
		}
		out2 = getFDR(PostProbSave)
			
		mylist = list(z=z, localfdr=out2$localfdr, FDR=out2$FDR, X=X, grid=l1$mids, breaks=l1$breaks,
        	grid.fz=l1$zdens, grid.f0z=grid.f0z, grid.zcounts=l1$zcounts,
        	dnull = M0, dmix=MTot, dalt = M1, mu0=mu0, sig0=sig0, p0=p0, ncomps=ncomps,
	        priorprob = PriorProbSave, postprob = PostProbSave, m1 = M1Save,
	        betasave = BetaSave,
	        weights = WeightsSave, means=MuSave, vars = VarSave
	    )
	}
	return(mylist);
}
