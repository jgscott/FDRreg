\name{FDRreg-package}
\alias{FDRreg-package}
\docType{package}
\title{
False discovery rate regression
}
\description{
False discovery rate regression.  Fits a model for multiple testing whereby the local false discovery rate may depend upon covariates via a logistic regression.
}
\details{
\tabular{ll}{
Package: \tab FDRreg\cr
Type: \tab Package\cr
Version: \tab 0.2\cr
Date: \tab 2015-04-18\cr
License: \tab GPL (>=3) \cr
}
The workhouse function is FDRreg(z, features, ...), where z is an observed vector of z statistics, and features is a matrix of covariates.  Do not add a column of ones to features to get an intercept term; the function does that for you, just like R's base lm() and glm() functions.
}
\author{
Author: James G. Scott, with contributions from Rob Kass and Jesse Windle.

Maintainer: James G. Scott  <james.scott@mccombs.utexas.edu>
}


\keyword{ False discovery rates }

\examples{

library(FDRreg)

# Simulated data
P = 2
N = 10000
betatrue = c(-3.5,rep(1/sqrt(P), P))
X = matrix(rnorm(N*P), N,P)
psi = crossprod(t(cbind(1,X)), betatrue)
wsuccess = 1/{1+exp(-psi)}

# Some theta's are signals, most are noise
gammatrue = rbinom(N,1,wsuccess)
table(gammatrue)

# Density of signals
thetatrue = rnorm(N,3,0.5)
thetatrue[gammatrue==0] = 0
z = rnorm(N, thetatrue, 1)
hist(z, 100, prob=TRUE, col='lightblue', border=NA)
curve(dnorm(x,0,1), add=TRUE, n=1001)

\dontrun{
# Fit the model
fdr1 <- FDRreg(z, X, nulltype='theoretical')

# Show the empirical-Bayes estimate of the mixture density
# and the findings at a specific FDR level
Q = 0.1
plotFDR(fdr1, Q=Q)

# Fitted model
fdr1$model

# Compare actual versus estimated prior probabilities of being a signal
plot(wsuccess, fdr1$priorprob, log='xy')
abline(0,1)

# Covariate effects
plot(X[,1], log(fdr1$priorprob/{1-fdr1$priorprob}), ylab='Logit of prior probability')
plot(X[,2], log(fdr1$priorprob/{1-fdr1$priorprob}), ylab='Logit of prior probability')

# Plot local FDR versus z statistics
plot(z, fdr1$localfdr, ylab='Local false-discovery rate')

# Extract findings at level FDR = Q
myfindings = which(fdr1$FDR <= Q)
hist(z[myfindings], breaks=50, col='lightblue', border='blue')
table(truth = gammatrue, guess = {fdr1$FDR <= Q})
}

}

\references{
James G. Scott, Ryan C. Kelly, Matthew A. Smith, Pengcheng Zhou, and Robert E. Kass (2015).  False discovery rate regression: application to neural synchrony detection in primary visual cortex.    Journal of the American Statistical Association, DOI: 10.1080/01621459.2014.990973. arXiv:1307.3495 [stat.ME].

Efron (2004). Large-scale simultaneous hypothesis testing: the choice of a null hypothesis. J. Amer. Statist. Assoc. 99, 96-104.

Efron (2005). Local false discovery rates. Preprint, Dept. of Statistics, Stanford University.

}


