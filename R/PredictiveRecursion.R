# April 11, 2014: migrated this functionality from a separate package over here

prfdr = function(z, mu0=0.0, sig0=1.0, control=list()) {
	
	mycontrol=list(gridsize = 500, npasses=10, decay=-0.67)
	mycontrol[(namc <- names(control))] <- control
	stopifnot(mycontrol$decay < -2/3, mycontrol$decay > -1, mycontrol$npasses >= 1)
	if(mycontrol$npasses < 5) warning("It is not recommended to use predictive recursion with < 5 passes through the data.\n")
	
	# Initial guess for alternative density and pi0
	x_grid = seq(min(z), max(z), length=mycontrol$gridsize)
	l1 = efron(z)
	nullprob = l1$p0
	# weights = pmax(0, 1 - l1$p0*dnorm(z, mu0, sig0)/l1$fz)
	# d1 = density(z, weights=weights/sum(weights), 1)
	# ispl1 = splines::interpSpline(d1$y ~ d1$x)
 	# theta_guess = predict(ispl1, x_grid)$y

    theta_guess = 0.05+seq(-1,1, length=mycontrol$gridsize)^2
	theta_guess = (1.0-nullprob)* theta_guess/trapezoid(x_grid, theta_guess)

	# We sweep through the data npasses times in random order
	# Set up this vector
	N = length(z)
	zz = rep(0, mycontrol$npasses*N)
	for(k in 1:mycontrol$npasses) {
		zz[1:N + (k-1)*N] = sample(z)
	}
	
	out1 = PredictiveRecursionFDR(zz, x_grid, theta_guess, nullprob=nullprob,
		mu0=mu0, sig0=sig0, decay=mycontrol$decay)
	out2 = eval_pr_dens(z, x_grid, out1$theta_subdens, sig0)
	# Use spline interpolation to get the mixture density estimates at the observed points
	ispl3 = splines::interpSpline(out1$y_mix ~ x_grid)
	fmix_z = predict(ispl3, z)$y
	
	out2 <- list(x_grid = x_grid, pi0 = out1$pi0, ftheta_grid = out1$theta_subdens, 
				fsignal_grid = out1$y_signal, fmix_grid = out1$y_mix,
				fsignal_z=out2$fsignal_z, fmix_z=fmix_z)
	out2
}
