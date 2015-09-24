# April 11, 2014: migrated this functionality from a separate package over here

prfdr = function(z, mu0=0.0, sig0=1.0, control=list()) {
	
	mycontrol=list(gridsize = 500, npasses=10, decay=-0.67)
	mycontrol[(namc <- names(control))] <- control
	stopifnot(mycontrol$decay < -2/3, mycontrol$decay > -1, mycontrol$npasses >= 1)
	if(mycontrol$npasses < 5) warning("It is not recommended to use predictive recursion with < 5 passes through the data.\n")
	
	# Initial guess for alternative density and pi0
	x_grid = seq(min(z), max(z), length=mycontrol$gridsize)
	nullprob = 0.95
	# theta_guess = pmin(1, 0.1+seq(-2,2, length=mycontrol$gridsize)^2)
	theta_guess = rep(1, length(x_grid))
	theta_guess = (1.0-nullprob)* theta_guess/trapezoid(x_grid, theta_guess)

	# We sweep through the data npasses times in random order
	# Set up the vector of indices
	N = length(z)
	sweeporder = rep(0L, mycontrol$npasses*N)
	for(k in 1:mycontrol$npasses) {
		sweeporder[1:N + (k-1)*N] = sample(0:(N-1))
	}
	
	out1 = PredictiveRecursionFDR(z, sweeporder, x_grid, theta_guess,
		mu0, sig0, nullprob=nullprob, decay=mycontrol$decay)
	out2 = eval_pr_dens(z, mu0, rep(sig0, N), x_grid, out1$theta_subdens)
	
	# pi(theta), the prior (normalized)
	pitheta_grid = out1$theta_subdens / (1-out1$pi0)

	# The mixture density
	f1_z = out2$fsignal_z
	f0_z = dnorm(z, mu0, sig0)
	fmix_z = (1.0-out1$pi0) * f1_z + out1$pi0*f0_z
	postprob = (1.0-out1$pi0) * f1_z / fmix_z
	
	out3 <- list(x_grid = x_grid, pi0 = out1$pi0, pitheta_grid = pitheta_grid, 
				f1_grid = out1$y_signal, fmix_grid = out1$y_mix,
				f0_z = f0_z, f1_z = out2$fsignal_z, fmix_z = fmix_z, postprob = postprob)
	out3
}


prfdr_het = function(z, mu0=0.0, sig0, control=list()) {
	# Predictive recursion estimate of a mixing density under heteroscedastic Gaussian error

	# Set up control parameters
	mycontrol=list(gridsize = 500, npasses=20, decay=-0.67)
	mycontrol[(namc <- names(control))] <- control
	stopifnot(mycontrol$decay < -.5, mycontrol$decay > -1, mycontrol$npasses >= 1)
	if(mycontrol$npasses < 5) warning("It is not recommended to use predictive recursion with < 5 passes through the data.\n")
	
	# Initial guess for alternative density and pi0
	x_grid = seq(min(z), max(z), length=mycontrol$gridsize)
	nullprob = 0.95
	# theta_guess = pmin(1, 0.1+seq(-2,2, length=mycontrol$gridsize)^2)
	theta_guess = rep(1, length(x_grid))
	theta_guess = (1.0-nullprob)* theta_guess/trapezoid(x_grid, theta_guess)

	# We sweep through the data npasses times in random order
	# Set up the vector of indices
	N = length(z)
	sweeporder = rep(0L, mycontrol$npasses*N)
	for(k in 1:mycontrol$npasses) {
		sweeporder[1:N + (k-1)*N] = sample(0:(N-1))
	}
	
	# Call cpp routine
	out1 = PredictiveRecursion_DifferentSigma(z, mu0, sig0, sweeporder, x_grid, theta_guess,
		nullprob=nullprob, decay=mycontrol$decay)
	out2 = eval_pr_dens(z, mu0, sig0, x_grid, out1$theta_subdens)
	
	out3 <- list(x_grid = x_grid, pi0 = out1$pi0, ftheta_grid = out1$theta_subdens, 
				fsignal_z=out2$fsignal_z)
	out3
}
