FDRsmooth1D = function(z, x=NULL, lambda=NULL,
	nulltype='theoretical', method='pr', stderr = NULL, control=list()) {
# False discovery rate smoothing
# z = vector of z scores
# x = integer locations of entries in z along 1D chain graph
# 		if x is missing, z is assumed to be in order, i.e. x=(1,...,N)
# nulltype = flag for what kind of null hypothesis to assume, theoretical/empirical/heteroscedastic
	
	stopifnot(method=='pr')
	
	# Set up control parameters
	mycontrol = list()
	if(method=='pr') {
		mycontrol$reltol = 1e-6
		mycontrol$maxit = 1000
		mycontrol$gridsize = 300
		mycontrol$decay = -0.67
		mycontrol$npasses = 10
		mycontrol$densknots=10
	}
	# Overwrite with user choices
	mycontrol[(namc <- names(control))] <- control
	
	# Sample size
	n = length(z)
	
	
	# Deal with ordering of z scores
	if(missing(x)) {
		z_sort = z
		x_order = 1:n
	} else {
		x_order = order(x)
		z_sort = z[x_order]
	}
	
	
	# Estimate the marginal density
	if(method=='pr') {
		
		# Compute f0 and f1, the marginals under null and alternative for each observation
		if(nulltype=='empirical') {

			# Currently using my implementation of Efron's central matching estimator
			l1 = efron(z, nmids=150, df=mycontrol$densknots, nulltype=nulltype)
			mu0 = l1$mu0
			sig0 = l1$sig0

			prfit = prfdr(z_sort, mu0, sig0, control=mycontrol)
			fmix_grid = prfit$fmix_grid
			f0_grid = dnorm(prfit$x_grid, mu0, sig0)
			f1_grid = prfit$f1_grid
		} else if(nulltype=='heteroscedastic') {
			if(missing(stderr)) {
				stop("Must specify standard error (stderr) if assuming heteroscedastic null.")
			}
			mu0 = 0.0
			sig0 = stderr
			prfit = prfdr_het(z_sort, mu0, sig0[x_order], control=mycontrol)
			fmix_grid = NULL
			f0_grid = NULL
			f1_grid = NULL
		} else {
			mu0 = 0.0
			sig0 = 1.0
			prfit = prfdr(z_sort, mu0, sig0, control=mycontrol)
			fmix_grid = prfit$fmix_grid
			f0_grid = dnorm(prfit$x_grid, mu0, sig0)
			f1_grid = prfit$f1_grid
		}
		# Extract marginal densities and fit regression
		p0 = prfit$pi0
		f0 = dnorm(z_sort, mu0, sig0)
		f1 = prfit$f1_z
		f1zeros = which(f1 < .Machine$double.eps)
		if(length(f1zeros > 0)) {
			f1[f1zeros] = min(f1[-m1zeros]) # substitute in the smallest nonzero value
			f1 = pmax(f1, .Machine$double.eps) # shouldn't happen but just in case!
		}
		x_grid = prfit$x_grid
	}

	# Now run FDR smoothing
	if(missing(lambda)) {
		lambda_grid = 10^seq(2,-0.5,length=50)
	} else {
		lambda_grid = lambda
	}
	n_lambda = length(lambda_grid)
	dof_grid = rep(0, n_lambda)
	dev_grid = rep(0, n_lambda)
	theta = matrix(0, nrow=n, ncol=n_lambda)
	theta_start = rep(0, n)
	fit_list = list()
	for(i in seq_along(lambda_grid)) {
		lambda = lambda_grid[i]
		this_fit = fdrs1D(f0, f1, lambda=lambda, mycontrol$reltol, mycontrol$maxit)
		dof_grid[i] = this_fit$n_jumps
		dev_grid[i] = 2*this_fit$val
		# update warm-start value and store result
		theta_start = this_fit$beta_hat
		theta[,i] = theta_start
	}
	bic_grid = dev_grid + sqrt(n) * dof_grid
	aic_grid = dev_grid + 2 * dof_grid
	list(theta = theta, prfit = prfit,
		lambda = lambda_grid, dof = dof_grid, dev = dev_grid, bic = bic_grid, aic = aic_grid)
}
	
	
fdrs1D = function(f0, f1, lambda, reltol, maxit, beta_start = NULL) {
# Used internally by the main wrapper function FDRsmooth1D
# f0 = vector of densities of z[i] under null
# f1 = vector of densities of z[i] under alternative
	
	n = length(f0)
	if(missing(beta_start)) {
		beta_hat = rep(0, n)
	} else {
		beta_hat = beta_start
	}
	
	prior_prob = ilogit(beta_hat)
	objective_old = -sum(log(prior_prob*f1 + (1-prior_prob)*f0))

	# Start iterates
	drift = max(1, 10*reltol)
	step_counter = 0
	while( {drift > reltol} & {step_counter <= maxit} ) {
		
		# E step
		prior_prob = ilogit(beta_hat)
		m1 = prior_prob*f1
		m0 = (1-prior_prob)*f0
		post_prob = m1/(m1+m0)
		
		# M step
		ebeta = exp(beta_hat)
		weights = ebeta/{(1+ebeta)^2}
		y = {(1+ebeta)^2}*post_prob/ebeta + beta_hat - (1+ebeta)
		weights = prior_prob*(1-prior_prob)
		y = beta_hat - (prior_prob - post_prob)/weights
		
		# Fit the 1-D fused lasso using the GLM pseudo-responses and weights
		beta_hat = fl_dp_weight(drop(y), w=drop(weights), lam=lambda)
		prior_prob = ilogit(beta_hat)
	
		objective_new = -sum(log(prior_prob*f1 + (1.0-prior_prob)*f0))
		drift = abs(objective_old - objective_new)/(abs(objective_new) + reltol)
		objective_old = objective_new
		step_counter = step_counter + 1
	}

	# Finish up and return fitted prior/posterior probabilities
	n_jumps = sum(abs(diff(beta_hat)) > 1e-6)
	m1 = prior_prob*f1
	m0 = (1-prior_prob)*f0
	post_prob = m1/(m1+m0)
	list(prior_prob = prior_prob, post_prob = post_prob,
		beta_hat = beta_hat, n_jumps = n_jumps, val = objective_old)
}

