\name{FDRreg-package}
\alias{FDRreg-package}
\docType{package}
\title{
False discovery rate regression
}
\description{
Tools for FDR problems, including false discovery rate regression.  Fits models whereby the local false discovery rate may depend upon covariates, either via a linear or additive logistic regression model.
}
\details{
\tabular{ll}{
Package: \tab FDRreg\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-02-25\cr
License: \tab GPL (>=3) \cr
}
The workhouse function is FDRreg(z,X, ...), where z is an observed vector of z statistics, and X is a matrix of covariates.  Do not add a column of ones to X to get an intercept term; the function does that for you, just like R's base lm() and glm() functions.
}
\author{
Author: James G. Scott, with contributions from Rob Kass and Jesse Windle.

Maintainer: James G. Scott  <james.scott@mccombs.utexas.edu>
}

\references{
False discovery rate regression: application to neural synchrony detection in primary visual cortex.  James G. Scott, Ryan C. Kelly, Matthew A. Smith, Pengcheng Zhou, and Robert E. Kass. arXiv:1307.3495 [stat.ME].
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
fdr1 <- FDRreg(z, covars=X, nmc=2500, nburn=100, nmids=120, nulltype='theoretical')


# Show the empirical-Bayes estimate of the mixture density
# and the findings at a specific FDR level
Q = 0.1
plotFDR(fdr1, Q=Q, showfz=TRUE)

# Posterior distribution of the intercept
hist(fdr1$betasave[,1], 20)

# Compare actual versus estimated prior probabilities of being a signal
plot(wsuccess, fdr1$priorprob)

# Covariate effects
plot(X[,1], log(fdr1$priorprob/{1-fdr1$priorprob}), ylab='Logit of prior probability')
plot(X[,2], log(fdr1$priorprob/{1-fdr1$priorprob}), ylab='Logit of prior probability')

# Local FDR
plot(z, fdr1$localfdr, ylab='Local false-discovery rate')

# Extract findings at level FDR = Q
myfindings = which(fdr1$FDR <= Q)
}

}

