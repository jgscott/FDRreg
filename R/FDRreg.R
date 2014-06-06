FDRreg = function(z, features, nulltype='theoretical', method='pr', stderr = NULL, control=list()) {
# False discovery rate regression
# z = vector of z scores
# X = design matrix of covariates, assumed NOT to have an intercept just as in vanilla lm()
# nulltype = flag for what kind of null hypothesis to assume, theoretical/empirical/heteroscedastic
	
	stopifnot(any(method=='pr', method=='efron'))
	
	# Set up control parameters
	mycontrol = list(center=TRUE, scale=TRUE)
	if(method=='pr') {
		mycontrol$gridsize = 300
		mycontrol$decay = -0.67
		mycontrol$npasses = 10
		mycontrol$lambda = 0.01
		mycontrol$densknots=10
	} else if(method=='efron') {
		mycontrol$gridsize = 150
		mycontrol$nmids=150
		mycontrol$densknots=10
	}
	# Overwrite with user choices
	mycontrol[(namc <- names(control))] <- control
	
	# Matrix of regressors, centered and scaled as requested
	N = length(z)
	X = scale(features, center=mycontrol$center, scale=mycontrol$scale)
	P = ncol(X)
	
	# Estimate the marginal density
	if(method=='pr') {
		
		# Compute M0 and M1, the marginals under null and alternative for each observation
		if(nulltype=='empirical') {

			# Currently using my implementation of Efron's central matching estimator
			l1 = efron(z, nmids=150, df=mycontrol$densknots, nulltype=nulltype)
			mu0 = l1$mu0
			sig0 = l1$sig0

			prfit = prfdr(z, mu0, sig0, control=mycontrol)
			fmix_grid = prfit$fmix_grid
			f0_grid = dnorm(prfit$x_grid, mu0, sig0)
			f1_grid = prfit$f1_grid
		} else if(nulltype=='heteroscedastic') {
			if(missing(stderr)) {
				stop("Must specify standard error (stderr) if assuming heteroscedastic null.")
			}
			mu0 = 0.0
			sig0 = stderr
			prfit = prfdr_het(z, mu0, sig0, control=mycontrol)
			fmix_grid = NULL
			f0_grid = NULL
			f1_grid = NULL
		} else {
			mu0 = 0.0
			sig0 = 1.0
			prfit = prfdr(z, mu0, sig0, control=mycontrol)
			fmix_grid = prfit$fmix_grid
			f0_grid = dnorm(prfit$x_grid, mu0, sig0)
			f1_grid = prfit$f1_grid
		}

		# Extract marginal densities and fit regression
		p0 = prfit$pi0
		M0 = dnorm(z, mu0, sig0)
		M1 = prfit$f1_z
		m1zeros = which(M1 < .Machine$double.eps)
		if(length(m1zeros > 0)) {
			M1[m1zeros] = min(M1[-m1zeros]) # substitute in the smallest nonzero value
			M1 = pmax(M1, .Machine$double.eps) # shouldn't happen but just in case!
		}
		x_grid = prfit$x_grid
		regressfit = fdrr_regress_pr(M0, M1, X, 1-p0, lambda=mycontrol$lambda)

	} else if(method=='efron') {
		if(nulltype=='heteroscedastic') {
			stop("Cannot use Efron's method under a heteroscedastic null.")
		}
		l1 = efron(z, nmids=mycontrol$nmids, df=mycontrol$densknots, nulltype=nulltype)
		mu0 = l1$mu0
		sig0 = l1$sig0
		p0 = l1$p0
		M0 = dnorm(z, mu0, sig0)
		M1 = NULL
		MTot = l1$fz	
		x_grid = l1$mids
		fmix_grid = l1$zdens
		f0_grid = dnorm(x_grid, mu0, sig0)
		f1_grid = NULL
		regressfit = fdrr_regress_efron(M0, MTot, X, 1-p0, N)
	}

	out2 = getFDR(regressfit$PostProb)
	list(	z=z, X=X, localfdr=out2$localfdr, FDR=out2$FDR, x_grid = x_grid,
			M0 = M0, M1 = M1,
			fmix_grid=fmix_grid, f0_grid = f0_grid, f1_grid = f1_grid, 
			mu0=mu0, sig0=sig0, p0=p0, priorprob = regressfit$W,
			postprob = regressfit$PostProb, model=regressfit$model
    )

}



BayesFDRreg = function(z, features, mu0=NULL, sig0 = NULL, empiricalnull=FALSE, nmc=5000, nburn=1000,
	control=list(), ncomps=NULL, priorpars = NULL) {
# Fully Bayesian version of false discovery rate regression
# z = vector of z scores
# features = design matrix of covariates, assumed NOT to have an intercept just as in vanilla lm()
# nulltype = flag for what kind of null hypothesis to assume, theoretical or empirical
# ncomps = how many mixture components for the alternative hypothesis

	mycontrol = list(center=FALSE, scale=FALSE, verbose=nmc+nburn+1)
	mycontrol[(namc <- names(control))] <- control
	
	N = length(z)
	X = cbind(1,scale(features, center= mycontrol$center, scale= mycontrol$scale))
	P = ncol(X)
	
	if(empiricalnull) {
		l1 = efron(z, nmids=150, df=15, nulltype='empirical')
		mu0 = l1$mu0
		sig0 = l1$sig0
		p0 = l1$p0
	} else {
		if(missing(sig0)) sig0 = rep(1,N)
		if(missing(mu0)) mu0 = 0
		p0 = NULL
	}
	sig0squared = sig0^2
	M0 = dnorm(z, mu0, sig0)

	# Initialize MCMC
	if(missing(priorpars)) {
		PriorPrec = diag(rep(1/25, P))
		PriorMean = rep(0,P)
	} else{
		PriorPrec = priorpars$PriorPrec
		PriorMean = priorpars$PriorMean
	}
		
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
		M1 = dnormix(z, emfit$weights[1:ncomps]/sum(emfit$weights[1:ncomps]),
						emfit$means[1:ncomps], emfit$vars[1:ncomps])
	} else 	M1 = dnorm(z, 0, 4)
	
	PriorPrecXMean = PriorPrec %*% PriorMean
	Beta = rep(0,P)
	Beta[1] = -3
	BetaSave = matrix(0, nrow=nmc, ncol=P)
	MuSave = matrix(0, nrow=nmc, ncol=ncomps)
	VarSave = matrix(0, nrow=nmc, ncol=ncomps)
	WeightsSave = matrix(0, nrow=nmc, ncol=ncomps)
	PostProbSave = 0
	PriorProbSave = 0
	M1Save = 0
	
	# Alternative hypothesis
	comp_weights = rep(1/ncomps, ncomps)
	comp_means = quantile(z[abs(z/sig0)>2], probs=seq(0.025,0.975,length=ncomps))
	comp_variance = rep(1, ncomps)
	myvar = comp_variance + 1
	
	# Main MCMC
	for(t in 1:(nmc+nburn)) {
		
		if(t %% mycontrol$verbose == 0) cat(t, "\n")
		
		### Update indicators
		Psi = drop(X %*% Beta)
		W = ilogit(Psi)
		PostProb = W*M1/{(1-W)*M0 + W*M1}		
		Gamma = rbinom(N,1,PostProb)
		
		
		### Update mixture of normals model
		
		# # Sufficient statistics
		# signals = z[which(Gamma==1)]
		# components = draw_mixture_component(signals, sig0, weights=comp_weights, mu = comp_means, tau2 = comp_variance) + 1
		# sumsignals = mosaic::maggregate(signals ~ factor(components, levels=1:ncomps), FUN='sum')
		# nsig = mosaic::maggregate(signals ~ factor(components, levels=1:ncomps), FUN='length')
		# tss = mosaic::maggregate((signals-comp_means[components])^2 ~ factor(components, levels=1:ncomps), FUN='sum')
		
		# #print(rbind(nsig, sumsignals, tss))

		# # Updates
		# for(k in 1:ncomps) myvar[k] = 1.0/rtgamma_once({nsig[k]+2}/2, {tss[k]+2}/2, 0, 1-.Machine$double.eps)
		# comp_variance = myvar - 1.0
		# muvar = comp_variance/{nsig + comp_variance*0.01}
		# muhat = sumsignals/{nsig + comp_variance*0.01}
		# comp_means = rnorm(ncomps, muhat, sqrt(muvar))	
		# comp_weights = rdirichlet_once(rep(1, ncomps) + nsig)

		cases = which(Gamma==1)
		signals = z[cases]
		components = draw_mixture_component(signals, sig0[cases], weights=comp_weights, mu = comp_means, tau2 = comp_variance) + 1

		# Draw latent means
		if(length(cases) > 0) {
			latentmeans.var = 1.0/(1.0/sig0[cases]^2 + 1.0/comp_variance[components])
			latentmeans.mu = latentmeans.var*(signals/(sig0[cases]^2) + (comp_means/comp_variance)[components])
			latentmeans = rnorm(length(cases), latentmeans.mu, sqrt(latentmeans.var))
			nsig = mosaic::maggregate(signals ~ factor(components, levels=1:ncomps), FUN='length')
			tss_thetai = mosaic::maggregate((latentmeans-comp_means[components])^2 ~ factor(components, levels=1:ncomps), FUN='sum')
			sum_thetai = mosaic::maggregate(latentmeans ~ factor(components, levels=1:ncomps), FUN='sum')
		} else {
			nsig = rep(0,ncomps)
			tss_thetai = rep(0,ncomps)
			mean_thetai = rep(0,ncomps)
		}

		# Updates
		for(k in 1:ncomps) comp_variance[k] = 1.0/rgamma(1, {nsig[k]+2}/2, rate={tss_thetai[k]+2}/2)
		muvar = comp_variance/{nsig + comp_variance*0.1}
		muhat = sum_thetai/{nsig + comp_variance*0.1}
		comp_means = rnorm(ncomps, muhat, sqrt(muvar))	
		comp_weights = rdirichlet_once(rep(5, ncomps) + nsig)
		M1 = marnormix(z, sig0squared, comp_weights, comp_means, comp_variance)

		### Update latent variables in logit likelihood
		Om = as.numeric(BayesLogit::rpg(N,rep(1,N),Psi))

		### Update regression parameters
		Kap = PostProb - 0.5
		PrecMat = t(X) %*% {Om * X} + PriorPrec
		Beta.V = solve(PrecMat)
		Beta.mu = Beta.V %*% {t(X) %*% Kap + PriorPrecXMean}
		Beta = t(mvtnorm::rmvnorm(1,mean=Beta.mu,sigma=Beta.V))	
		if(t > nburn) {
			BetaSave[t-nburn,] = Beta
			MuSave[t-nburn,] = comp_means
			VarSave[t-nburn,] = comp_variance
			WeightsSave[t-nburn,] = comp_weights
			PostProbSave = PostProbSave + (1.0/nmc)*PostProb
			PriorProbSave = PriorProbSave + (1.0/nmc)*W
			M1Save = M1Save + (1.0/nmc)*M1
		}
	}
	out2 = getFDR(PostProbSave)
		
	mylist = list(z=z, localfdr=out2$localfdr, FDR=out2$FDR, X=X,
			fnull_z = M0, fsignals_z = M1Save, mu0=mu0, sig0=sig0, p0=p0, ncomps=ncomps,
			priorprob = PriorProbSave, postprob = PostProbSave, 
			coefficients = BetaSave, weights = WeightsSave, means=MuSave, vars = VarSave
    )
	return(mylist);
}
