ilogit = function(x) 1/{1+exp(-x)}
flogit = function(x) log(x/{1-x})

# Truncated gamma random draws
rtgamma = function(n, a, b, lb, ub) {
	lub = pgamma(lb, a, rate=b)
	uub = pgamma(ub, a, rate=b)
	u = runif(n, lub, uub)
	qgamma(u, a, rate=b)
}


SoftLogitFit = function(y, X, lambda=1e-6, start=NULL) {
	# y: a vector of "fuzzy" ones and zeros, i.e. expected binary outcomes
	# X: design matrix assumed to include a column of 1's for an intercept
	# lambda is a ridge penalty parameter, defaulting to 1e-6
	#	(basically no regularization, only there to ensure numerical stability)
	if(missing(start)) {
		start = rep(0, ncol(X))
	}
	mymax = optim(start, fn = SoftLogitLoss, gr = SoftLogitGradient,
		method='BFGS', y=y, X = X, lambda=lambda, hessian=TRUE)
	list(coef = mymax$par, value = mymax$value, hessian = mymax$hessian)
}


efron = function(z, nmids=150, pct=-0.01, pct0=0.25, df=10, nulltype='theoretical') {
	# estimate f(z) and f_0(z) using Efron (2004)'s method
	stopifnot(any(nulltype == 'theoretical', nulltype=='empirical'))
	
	result = tryCatch({

		N = length(z)
		med = median(z)
		myrange = med + (1 - pct) * (range(z) - med)
		lb = myrange[1]
		ub = myrange[2]
		                
		breaks = seq(lb, ub, length= nmids +1)
		h1 = hist(z, breaks = breaks, plot = FALSE)
		mids = (breaks[-1] + breaks[-(nmids+1)])/2
		zcounts = h1$counts
		glm1 = glm(zcounts ~ splines::ns(mids, df = df), family=poisson)
		zrate = glm1$fit
		D = (zcounts - zrate)/sqrt(zrate+1)
		D = D[-c(1,nmids)]
		if (sum(D^2) > qchisq(0.9, nmids-2-df)) {
			warning(paste0("f(z) misfit = ", round(D, 1), ".  Rerun with increased df."))
		}	
		
		zdens = {zrate/sum(zrate)}/diff(breaks)

		# Now do spline interpolation for the density at the observed points
		ispl2 = splines::interpSpline( zdens ~ mids )
		fz = predict(ispl2, z)$y
		
		# Pick out the middle of the data points
		ql = quantile(z, pct0)
		qu = quantile(z, 1-pct0)
		ind0 = intersect(which(z > ql), which(z<qu))
		if(nulltype=='empirical') {
		# empirical null by central moment matching
			z0 = z[ind0]
			l0 = log(fz[ind0])
			zmax = z[which.max(l0)]
			lm0 = lm(l0~I(z0-zmax) + I((z0-zmax)^2))
			b0 = coef(lm0)
			sig = as.numeric(sqrt(-1.0/{2.0*b0[3]}))
			mu = as.numeric(-b0[2]/(2.0*b0[3]) + zmax)
			# lm0 = lm(l0 ~ z0 + I(z0^2))
			# b0 = coef(lm0)
			# sig = as.numeric(sqrt(-1/{2*b0[3]}))
			# mu = as.numeric(b0[2] * sig^2)
		} else {
		# theoretical null
			sig = 1
			mu = 0
		}
		p0 = sum(fz[ind0])/sum(dnorm(z[ind0], mu, sig))
		localfdr = pmin(1, p0*dnorm(z, mu, sig)/fz)
		list(mids=mids, breaks=breaks, zcounts=zcounts, zdens=zdens,
			z=z, fz=fz, mu0=mu, sig0=sig, p0=p0, fdr=localfdr)
	}, error = function(err) {
		print(err)
		list(mids=NULL, breaks=NULL, zcounts=NULL, zdens=NULL,
			z=NULL, fz=NULL, mu0=NULL, sig0=NULL, p0=NULL, fdr=NULL)
	}, finally = {
		# Nothing to clean up
	})
	return(result)
}

getFDR = function(postprob) {
	# postprob is a vector of posterior probabilities
	# from which local fdr and (Bayesian) FDR are extracted
	indices = 1:length(postprob)
	iorder = order(postprob, decreasing=TRUE)
	porder = postprob[iorder]
	localfdr.order = 1-porder
	FDR.order = cumsum(localfdr.order)/indices
	localfdr = indices  # placeholder
	localfdr[iorder] = localfdr.order
	FDR = indices  # placeholder
	FDR[iorder] = FDR.order
	
	# Where local fdr is 1, report the most conservative FDR 
	fdrmax = which(localfdr == 1)
	FDR[fdrmax] = max(FDR)
	list(localfdr=localfdr, FDR=FDR)
}

plotFDR = function(fdrr, Q=0.1, showrug=TRUE, showsub=TRUE, breaks=150, ...) {
	N = length(fdrr$z)
	mytitle = paste0('')
	par(mar=c(5,4,1,1))
	par(...)
	h1 = hist(fdrr$z, breaks, plot=FALSE)
	plot(h1, freq=FALSE, col='lightgrey', border='grey',
		main=mytitle, xlab='', ylab='', axes=FALSE, ylim=c(0, 1.05*max(h1$density)))
	mysub = paste0('Grey bars: original z scores\nBlue bars: fraction signals in each bin')
	axis(1, pos=0, tick=FALSE, cex.axis=0.9)
	axis(2, tick=FALSE, las=1, cex.axis=0.9)
	zcut = data.frame(prob=fdrr$postprob, bucket=cut(fdrr$z, h1$breaks))
	pmean = mosaic::maggregate(prob~bucket, data=zcut, FUN=mean)
	pmean[is.na(pmean)] = 0
	par(new=TRUE)
	h2 = h1
	h2$density = h2$density * pmean
	plot(h2, freq=FALSE, col='blue', border='grey', axes=FALSE, xlab='', ylab='', main='', ylim=c(0, 1.05*max(h1$density)))
	lines(fdrr$x_grid, fdrr$fmix_grid, col='black')
	lines(fdrr$x_grid, fdrr$p0 * fdrr$f0_grid, col='red', lty='dotted', lwd=1)
	legend('topright', c(expression(f(z)), expression(pi[0] %.% f[0](z))),
		lty=c('solid', 'dotted'), col=c('black', 'red'), bty='n')
	if(showrug) {
		rug( fdrr$z[fdrr$FDR < Q], ticksize=0.03, col='black')
		mysub = paste0(mysub, '\nBlack rug: discoveries at FDR = ', Q)
	}
	if(showsub) title(sub=mysub)
	par(new=FALSE)
}


# Benjamini-Hochberg with two-sided p-values
BenjaminiHochberg = function(zscores, fdr_level) {
# zscores is a vector of z scores
# fdr_level is the desired level (e.g. 0.1) of control over FDR
# returns a binary vector where 0=nofinding, 1=finding at given FDR level
	N = length(zscores)
	pval2 = 2*pmin(pnorm(zscores), 1- pnorm(zscores))
	cuts = (1:N)*fdr_level/N
	bhdiff = sort(pval2)-cuts
	bhcutind2 = max(which(bhdiff < 0))
	bhcut2 = sort(pval2)[bhcutind2]
	0+{pval2 <= bhcut2}
}

# Utility function for extracting error rates and a confusion matrix
GetErrorRates = function(truth, guess) {
# truth is a binary vector saying which cases are signals
# guess is a binary vector saying which cases are "findings" from a given procedure
	confusion_matrix = table(factor(truth, levels=c(0,1)), factor(guess, levels=c(0,1)))
	
	# Need to catch the corner case: no signals
	if( sum(truth) == 0 ) {
		true_positive_rate = 0
	} else {
		true_positive_rate = confusion_matrix[2,2]/sum(truth)
	}

	# Need to catch the corner case: no discoveries
	if( sum(guess) == 0 )  {
		false_discovery_rate = 0
	} else {
		false_discovery_rate = confusion_matrix[1,2]/sum(guess)
	}

	list(tpr = true_positive_rate, fdr = false_discovery_rate, confusion = confusion_matrix)
}


# Iteratively fit the regression piece of the PR-based FDR regression
fdrr_regress_pr = function(M0, M1, X, W_initial, maxit=2500, abstol = 1e-6, lambda=0.01) {
	stopifnot( M0 > 0, M1 > 0 )
	P = ncol(X)
	result = tryCatch({

		# Initialize the regression fit
		travel=1
		Xs = cbind(1,X)
		PostProb = W_initial*M1/(W_initial*M1 + (1-W_initial)*M0)
		suppressWarnings(lm1 <- SoftLogitFit(PostProb, Xs, lambda=lambda))
		betaguess = lm1$coef
		W = ilogit(drop(Xs %*% betaguess))
		PostProb = W*M1/(W*M1 + (1-W)*M0)
		oldval = lm1$value

		# Iterate until convergence, each time with a warm start
		passcounter = 0
		while(abs(travel/oldval) > abstol && passcounter <= maxit) {
			suppressWarnings(lm1 <- SoftLogitFit(PostProb, Xs, lambda=lambda, start=betaguess))
			newval = lm1$value
			travel = abs((oldval - newval)/(oldval + abstol))
			oldval = newval
			betaguess = lm1$coef
			W = ilogit(drop(Xs %*% betaguess))
			PostProb = W*M1/(W*M1 + (1-W)*M0)
			passcounter = passcounter + 1
		}
		if(abs(travel/oldval) > abstol) {
			mywarning = paste0('\nMaximum FDRR iteration (maxit) reached for PR method. ',
						'Try re-running with a weaker tolerance or larger maxit.')
			warning(mywarning, immediate.=FALSE)
		}
		list(PostProb = PostProb, W = W, model = lm1)
	}, error = function(err) {
		print(err)
		print("An error was encountered in fitting the regression.  Reverting to the PR no-covariates model.")
		list(PostProb = W_initial*M1/(W_initial*M1 + (1-W_initial)*M0), W = W_initial, model = NULL)
	}, finally = {
		# No cleanup necessary
	})
	return(result)
}

# Iteratively fit the regression piece of the FDR regression using Efron's estimator
fdrr_regress_efron = function(M0, MTot, X, W_initial, maxit=500, abstol = 1e-6) {
	stopifnot( M0 > 0, MTot > 0 )
	P = ncol(X)
	result = tryCatch({

		# Initialize the regression fit
		travel=1
		nullcases = which((1-W_initial)*M0>MTot)	
		PostProb = pmax(0, pmin(1, 1 - (1-W_initial)*M0/MTot))
		PostProb[nullcases] = 0
		suppressWarnings(lm1 <- glm(PostProb ~ X, family=binomial))
		W = fitted(lm1)
		PostProb = pmax(0, pmin(1, 1 - (1-W)*M0/MTot))
		PostProb[nullcases] = 0
		betaguess = coef(lm1)


		# Iteratively refine estimate for beta
		travel=1
		passcounter = 0
		while(travel > abstol && passcounter <= maxit) {
			suppressWarnings(lm1 <- glm(PostProb ~ X, family=binomial, start=betaguess))
			newbeta = coef(lm1)
			W = fitted(lm1)
			PostProb = pmax(0, pmin(1, 1 - (1-W)*M0/MTot))
			PostProb[nullcases] = 0
			travel = sum(abs(betaguess-newbeta))
			betaguess = newbeta
			passcounter = passcounter + 1
		}

		if(travel > abstol) {
			myerror = paste0('Maximum FDRR iteration (maxit) reached for Efron\'s method. ',
						'Try re-running with a weaker tolerance or larger maxit.')
			warning(mywarning, immediate. = TRUE)
		}

		list(PostProb = PostProb, W = W, model = lm1)

	}, error = function(err) {
		#print(err)
		print("An error was encountered in fitting the regression.  Reverting to the Efron no-covariates model.")
		list(PostProb = pmax(0, pmin(1, 1 - (1-W_initial)*M0/MTot)), W = W_initial, model = NULL)
	}, finally = {
		# No cleanup necessary
	})
	return(result)
}




# EM deconvolution for a Gaussian mixture model for mu_i with N(mu_i,1) observations
deconvolveEM = function(z, ncomps, rel.tol=1e-7, plotit=FALSE) {
	N = length(z)
	
	comp_weights = c(0.5*rep(1/ncomps, ncomps), 0.5)
	comp_means = c(quantile(z[abs(z)>2], probs=seq(.025, 0.975, length=ncomps)), 0)
	comp_variance = c(rep(2, ncomps), 1)
	
	loglike = sum(log(dnormix(z, comp_weights, comp_means, comp_variance)))
	converged=FALSE
	
	marg_like = matrix(0, nrow=N, ncol=ncomps+1)
	marg_like[,ncomps+1] = dnorm(z)
	while(!converged) {
	
		# E step
		for(j in 1:ncomps) {
			marg_like[,j] = dnorm(z, comp_means[j], sqrt(comp_variance[j]))
		}
		post_probs = scale(marg_like, center=FALSE, scale=1/comp_weights)
		post_probs = post_probs/rowSums(post_probs)
		
		# M step
		comp_weights = colSums(post_probs)/N
		for(j in 1:ncomps) {
			comp_means[j] = sum(z*post_probs[,j])/sum(post_probs[,j])
			# The component-level variances are constrained to be >=1, consistent with a deconvolution
			comp_variance[j] = max(1,  sum( {(z-comp_means[j])^2} * post_probs[,j])/sum(post_probs[,j]) )
		}
		
		loglikenew = sum(log(dnormix(z, comp_weights, comp_means, comp_variance)))
		relative_change = abs(loglikenew - loglike)/abs(loglike + rel.tol)
		converged = {relative_change < rel.tol}
		loglike = loglikenew
	}
	if(plotit) {
		par(mfrow=c(1,2))
		hist(z, 100, prob=TRUE, col='grey', border=NA)
		curve(marnormix(x, rep(1, N), comp_weights, comp_means, comp_variance-1), add=TRUE, n=1000)
		hist(z, 100, prob=TRUE, col='grey', border=NA)
		curve(dnormix(x, comp_weights[1:ncomps]/sum(comp_weights[1:ncomps]), comp_means[1:ncomps], comp_variance[1:ncomps]), add=TRUE, n=1000)

	}
	list(weights = comp_weights, means = comp_means, vars = comp_variance, loglike = loglike, AIC = -2*loglike + 2*{3*ncomps -1 })
}
