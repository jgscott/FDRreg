FDRreg = function(z, covars, nulltype='empirical', type='linear', nmc=10000, nburn=500, nmids=150, densknots=10, regknots=5) {
# False discovery rate regression
# z = vector of z scores
# covars = design matrix of covariates, assumed NOT to have an intercept just as in vanilla lm()
# nulltype = flag for what kind of null hypothesis to assume, theoretical or empirical
	N = length(z)
	fjindex=list()
	if(type=='linear') {
		X = cbind(1, scale(covars, center=TRUE, scale=TRUE))
		fjindex = as.list(2:ncol(X))
	}
	else if(type=='additive') {
		X = matrix(rep(1,N))
		for(j in 1:ncol(covars)) {
			Xj = make.bspline.matrix(as.matrix(covars[,j]), nknots=regknots, method='equi')
			fjindex[[j]] = ncol(X) + 1:ncol(Xj)
			X = cbind(X, Xj)
		}
	}
	else {
		warning('Invalid model type specified: assuming linearity in the covariates.')
	}
	P = ncol(X)

	# Estimate the marginal density and empirical null
	l1 = efron(z, nmids=nmids, df=densknots)
	mu0 = l1$mu0
	sig0 = l1$sig0
	if(nulltype=='empirical') {
		M0 = dnorm(z,mu0,sig0)
		grid.f0z = dnorm(l1$mids,mu0,sig0)
	}
	else if(nulltype=='theoretical') {
		M0 = dnorm(z)
		grid.f0z = dnorm(l1$mids)
	}
	else {
		warning('Invalid null type specified: assuming empirical null.')
		M0 = dnorm(z,mu0,sig0)
		grid.f0z = dnorm(l1$mids,mu0,sig0)
	}
	MTot = l1$fz
	nullcases = which(M0>MTot)

	# Initialize MCMC
	PriorPrec = diag(c(1e-9, rep(0.1, P-1)))
	PriorMean = rep(0,P)
	
	# Pass to C++ and return result
	out1 = FDRregCPP(z, X, M0=M0, MTot = MTot, PriorPrec, PriorMean, nmc=nmc, nburn=nburn)
	out2 = getFDR(out1$postprob)
	list(z=z, localfdr=out2$localfdr, FDR=out2$FDR, X=X, grid=l1$mids, breaks=l1$breaks,
         grid.fz=l1$zdens, grid.f0z=grid.f0z, grid.zcounts=l1$zcounts, dnull = M0,
         dmix=MTot, empirical.null=list(mu0=mu0, sig0=sig0), betasave = out1$betasave,
         priorprob = out1$priorprob, postprob = out1$postprob, fjindex=fjindex)
}






BayesFDRreg = function(z, covars, nulltype='empirical', type='linear', nmc=10000, nburn=500, nmids=150, densknots=10, ncomps=5, regknots=5, verbose=500) {
# More Bayesian version of false discovery rate regression
# z = vector of z scores
# covars = design matrix of covariates, assumed NOT to have an intercept just as in vanilla lm()
# nulltype = flag for what kind of null hypothesis to assume, theoretical or empirical
	N = length(z)
	fjindex=list()
	if(type=='linear') {
		X = cbind(1, scale(covars, center=TRUE, scale=TRUE))
		fjindex = as.list(2:ncol(X))
	}
	else if(type=='additive') {
		X = matrix(rep(1,N))
		for(j in 1:ncol(covars)) {
			Xj = make.bspline.matrix(as.matrix(covars[,j]), nknots=regknots, method='equi')
			fjindex[[j]] = ncol(X) + 1:ncol(Xj)
			X = cbind(X, Xj)
		}
	}
	else {
		warning('Invalid model type specified: assuming linearity in the covariates.')
	}
	P = ncol(X)

	# Estimate the marginal density and empirical null
	l1 = efron(z, nmids=nmids, df=densknots)
	mu0 = l1$mu0
	sig0 = l1$sig0
	if(nulltype=='empirical') {
		M0 = dnorm(z,mu0,sig0)
		grid.f0z = dnorm(l1$mids,mu0,sig0)
	}
	else if(nulltype=='theoretical') {
		M0 = dnorm(z)
		grid.f0z = dnorm(l1$mids)
	}
	else {
		warning('Invalid null type specified: assuming empirical null.')
		M0 = dnorm(z,mu0,sig0)
		grid.f0z = dnorm(l1$mids,mu0,sig0)
	}
	MTot = l1$fz
	
	# Initialize MCMC
	PriorPrec = diag(c(1e-9, rep(0.1, P-1)))
	PriorMean = rep(0,P)
	
	# Pass to C++ and return result
	out1 = BayesFDRregCPP(z, X, M0=M0, ncomps, PriorPrec, PriorMean, nmc=nmc, nburn=nburn, grid=l1$mids)
	out2 = getFDR(out1$postprob)
	list(z=z, localfdr=out2$localfdr, FDR=out2$FDR, X=X, grid=l1$mids, breaks=l1$breaks,
         grid.fz=l1$zdens, grid.f0z=grid.f0z, grid.zcounts=l1$zcounts, dnull = M0, falt = out1$fthetasave,
         dmix=MTot, empirical.null=list(mu0=mu0, sig0=sig0), betasave = out1$betasave,
         musave = out1$musave, weightssave = out1$weightssave, varsave = out1$varsave,
         priorprob = out1$priorprob, postprob = out1$postprob, fjindex=fjindex)
         
	# # Set up the grid for the approximation to f1(z)
	# mids = l1$mids
	# breaks = l1$breaks

	# # Initialize MCMC
	# PriorPrec = diag(c(1e-9, rep(0.1, P-1)))
	# PriorPrecXMean = rep(0,P)
	# Beta = rep(0, P); Beta[1] = flogit(sum(z>2)/N)
	# BetaSave = matrix(0, nrow=nmc, ncol=P)
	# PostProbSave = 0
	# PriorProbSave = 0
	# M1Save = 0
	# #M1 = dnorm(z, mean(z[which(z>2)]), 2)
	# M1 = dnorm(z, 0, 4)
	# myvar = 16
	
	# # Alternative hypothesis
	# comp_weights = rep(1/ncomps, ncomps)
	# comp_means = quantile(z, probs=seq(1/ncomps, 1-1/ncomps, length=ncomps))
	# comp_variance = rep(1, ncomps)
	
	# # Main MCMC
	# for(t in 1:(nmc+nburn)) {
		
		# if(t %% verbose == 0) cat(t, "\n")
		
		# ### Update indicators
		# Psi = drop(X %*% Beta)
		# W = ilogit(Psi)
		# PostProb = W*M1/{(1-W)*M0 + W*M1}		
		# Gamma = rbinom(N,1,PostProb)
		
		# ### Update mixture of normals model
		
		# # Sufficient statistics
		# signals = z[which(Gamma==1)]
		# components = draw_mixture_component(signals, rep(1,N), weights=comp_weights, mu = comp_means, tau2 = comp_variance) + 1
		# sumsignals = mosaic::maggregate(signals ~ factor(components), FUN='sum')
		# nsig = mosaic::maggregate(signals ~ factor(components), FUN='length')
		# tss = mosaic::maggregate((signals-comp_means[components])^2 ~ factor(components), FUN='sum')
		
		# # Updates
		# myvar = 1/rtgamma(ncomps, {nsig}/2, {tss}/2, 0, 1)
		# comp_variance = myvar - 1
		# muvar = 1/{nsig/comp_variance + 0.01}
		# muhat = muvar*{sumsignals/comp_variance}
		# comp_means = rnorm(ncomps, muhat, sqrt(muvar))	
		# M1 = marnormix(z, rep(1,N), comp_weights, comp_means, comp_variance)
		

		# ### Update latent variables in logit likelihood
		# Om = as.numeric(rpg(N,rep(1,N),Psi))

		# ### Update regression parameters
		# Kap = PostProb - 1/2
		# PrecMat = t(X) %*% {Om * X} + PriorPrec
		# Beta.V = solve(PrecMat)
		# Beta.mu = Beta.V %*% {t(X) %*% Kap + PriorPrecXMean}
		# Beta = t(mvtnorm::rmvnorm(1,mean=Beta.mu,sigma=Beta.V))	
		# if(t > nburn) {
			# BetaSave[t-nburn,] = Beta
			# PostProbSave = PostProbSave + (1/nmc)*PostProb
			# PriorProbSave = PriorProbSave + (1/nmc)*W
			# M1Save = M1Save + (1/nmc)*M1
		# }
	# }
	# list(fz=M1Save, mu0=mu0, sig0=sig0,
		# X=X, BetaSave = BetaSave, PriorProb = PriorProbSave, PostProb = PostProbSave, fjindex=fjindex)
}
