########
# A simple example of Bayesian FDR using prfdr and prfdr_het
########
library(FDRreg)

# simulate from a mixture distribution of true means
N = 10000
w_true = 0.2
gamma_true = rbinom(N, 1, w_true)

# What is the alternative hypothesis?
sim_alt = function(n) {
	# unimodal
	rnorm(n, 0, 2)  
	
	# # bimodal
	# n1 = floor(n/2)
	# c(rnorm(n1, 3, 1), rnorm(n-n1, -3, 1))
}

theta_true = rep(0, N)
theta_true[gamma_true == 1] = sim_alt(sum(gamma_true))
plot(ecdf(theta_true))

# What are the standard errors?
sigma0 = sqrt(1/rgamma(N, 10, 10))
hist(sigma0)

# Simulate data and convert to z scores
y = rnorm(N, theta_true, sigma0)
z = y/sigma0

# Deconvolution under heteroskedastic model via predictive recursion
out1 = prfdr_het(y, 0, sigma0, control=list(npasses=10))

# Compare fit with truth
out1$pi0; 1-w_true
plot(out1$x_grid, out1$ftheta_grid/(1-out1$pi0), type='l')
lines(density(theta_true[gamma_true==1]), col='blue')

# Posterior probability of signal
w0 = out1$pi0*dnorm(y, 0, sigma0)
w1 = (1-out1$pi0)*out1$fsignal_z
postprob1 = w1/{w1+w0}

# Extract realized FDR
fdr_level = 0.1
FDR1 = getFDR(postprob1)
gamma1 = 0 + {FDR1$FDR <= fdr_level}

# Calculate FDR and true positive rates
GetErrorRates(gamma_true, gamma1)

# Compare z scores with posterior probabilities
plot(z, postprob1)

#####
# Compare with simple empirical-Bayes estimate under Scott and Berger (2006) prior
#####

# objective function for empirical Bayes by marginal MLE
target = function(parvec, y, sigma) {
	w = 1/{1+exp(-parvec[1])}
	tau2 = exp(parvec[2])
	loglike = sum(log( (1-w)*dnorm(y, 0, sigma) + w*dnorm(y, 0, sqrt(tau2+sigma^2)) ))
	-loglike
}


# Should not do especially well if truth is bimodal
# But should estimate the unimodal case OK
out3 = optim(c(0,0), target, y=y, sigma=sigma0)
w_hat = 1/{1+exp(-out3$par[1])}
tau = sqrt(exp(out3$par[2]))

# Calculate posterior probabilities
w0 = (1-w_hat)*dnorm(y, 0, sigma0)
w1 = w_hat*dnorm(y, 0, sqrt(tau^2 + sigma0^2))
postprob2 = w1/{w1+w0}
FDR2 = getFDR(postprob2)
gamma2 = 0 + {FDR2$FDR <= fdr_level}
GetErrorRates(gamma_true, gamma2)



#####
# Deconvolution of z scores instead
#####

mu_true = theta_true/sigma0
out2 = prfdr(z)
out2$pi0; 1-w_true

# Compare fit with truth
plot(out2$x_grid, out2$pitheta_grid, type='l')
lines(density(mu_true[gamma_true==1]), col='blue')


