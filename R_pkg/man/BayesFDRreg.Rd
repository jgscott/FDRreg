% package: FDRreg
\name{BayesFDRreg}

\alias{BayesFDRreg}

\title{Fully Bayesian False Discovery Rate Regression}

\description{
  Estimate a fully Bayes false-discovery rate regression model for test statistics z and regressors X.
}


\usage{
BayesFDRreg(z, features, mu0=NULL, sig0 = NULL, empiricalnull=FALSE, nmc=5000, nburn=1000,
  control=list(), ncomps=NULL, priorpars = NULL)
}

\arguments{
  \item{z}{An N dimensional vector; z_i is the test statistic for observation i.}
  \item{features}{An N x P dimensional design matrix; feature vector x_i for the ith observation is the ith row.  This is assumed not to have a column of ones representing an intercept.  Just like in R's own lm() and glm(), the intercept will be added by the fitting algorithm.  By default the features will be centered and scaled.  This can be modified using the 'control' flag.}
  \item{mu0}{The mean of the null density used to fit the model.  Defaults to NULL, in which case the 'empiricalnull' flag governs the choice of null hypothesis.}
  \item{sig0}{The standard deviation of the null density used to fit the model.  Can also be a vector of length N if each z[i] has a known standard error.  Defaults to NULL, in which case the 'empiricalnull' flag governs the choice of null hypothesis.}
  \item{empiricalnull}{If true, an empirical null is fitted using Efron's central matching estimator.  See Efron (2004).  If false and mu0/sig0 are specified, then these quantities are used for the null hypothesis.  If false and mu0/sig0 are missing, a standard normal null is assumed.}
  \item{nmc}{The number of MCMC iterations saved.  Defaults to 5000.}
  \item{nburn}{The number of initial MCMC iterations discarded as burn-in.  Defaults to 1000.}
  \item{ncomps}{The number of mixture components used in the Gaussian mixture model of the alternative density f^1(z).  Defaults to NULL; see Details.}
  \item{priorpars}{A list with elements PriorPrec (the prior precision matrix) and PriorMean (the prior mean vector).  If missing, the precision precision matrix is set 1/25 for every parameter, and the prior mean to 0.}
  \item{control}{A list of further options to the fitting algorithm; see Details below.}
}


\value{
  \item{z}{The test statistics provided as the argument z.}
  \item{X}{The design matrix used in the regression.  This will include an added column for an intercept and will reflect the requested centering and scaling of the original feature matrix.}
  \item{localfdr}{The vector of local false discovery rates (lfdr) corresponding to each element of z.  localfdr[i] is (1-p[i]). where p[i] is the fitted posterior probability that z[i] comes from the non-null (signal) population. Note localfdr is not necessarily monotonic in z, because the regression model allows the prior probability that z[i] is a signal to change with covariates X[i,].}
  \item{FDR}{The corresponding vector of cut-level false discovery rates (FDR) for the elements of z. Used for extracting findings at a given FDR level.  FDR[i] is the estimated false discovery rate for the cohort of test statistics whose local fdr's are at least as small as localfdr[i] --- that is, the z[j]'s such that localfdr[j] <= localfdr[i].  If you want the findings that meet a certain FDR cutoff, use this return value.}
  \item{M0}{The estimated (or assumed) null density at each of the observed z scores; M0[i] corresponds to f^0(z[i]).}
  \item{M1}{The estimated alternative density at each of the observed z scores; M1[i] corresponds to f^1(z[i]).}
  \item{mu0}{The mean of the null density used to fit the model.}
  \item{sig0}{The standard deviation of the null density used to fit the model.}
  \item{p0}{The estimated global null probability estimated by Efron (2004)'s method.  This is useful for comparison with the FDR regression results, because it assumes that the prior probability of being a signal does not change with covariates.}
  \item{ncomps}{The number of Gaussian mixture components used to fit the alternative hypothesis.}
  
 \item{priorprob}{The estimated prior probability of being a signa for each observation z_i.  Here priorprob[i] = P(z_i is non-null). }
 \item{postprob}{The estimated posterior probabilities of being a signal each observation z_i: postprob[i] = P(z_i is non-null | data), and localfdr[i] = 1-postprob[i]. }


 \item{coefficients}{Each row is a posterior draws of the regression vector. }
 \item{weights}{Each row is a posterior draw of the weights on the mixture components used to estimate the alternative hypothesis.}
 \item{means}{Each row is a posterior draw of the means of the mixture components used to estimate the alternative hypothesis.}
 \item{vars}{Each row is a posterior draw of the variances of the mixture components used to estimate the alternative hypothesis.}

}

\details{

  This model assumes that a z-statistic z arises from

  \deqn{ f(z_i) = w_i f^1(z) + (1-w_i) f^0(z) , }

where f^1(z) and f^0(z) are the densities/marginal likelihoods under the alternative and null hypotheses, respectively, and where w_i is the prior probability that z_i is a signal (non-null case).  Efron (2004)'s central-matching method may be used to estimate f^0(z) nonparametrically.  The prior probabilities w_i are estimated via logistic regression against covariates, using the Polya-Gamma Gibbs sampler of Polson, Scott, and Windle (JASA, 2013).  The alternative f^1(z) is estimated using a finite Gaussian mixture; see Details.


If ncomps=NULL (default), then ncomps will be chosen by sequentially adding a mixture component, fitting a two-groups model by EM assuming no covariate effects, and terminating when there is no further improvement to AIC.  This will not fix the model itself, merely the number of components for f^1(z) used by the FDR regression model, which is subsequently fit by MCMC.

  The control parameter can be a list with any of the following elements:
\itemize{
    \item center, scale: Should the feature matrix be centered and scaled?  Both default to true.
    \item verbose: How often (measured in MCMC iterations) should MCMC progress be reported?  Defaults to no progress output?
}

}

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
bfdr1 <- BayesFDRreg(z, X, empiricalnull=FALSE, nmc=2000, nburn=200,
        control=list(verbose=200))

# Show the empirical-Bayes estimate of the mixture density
# and the findings at a specific FDR level
Q = 0.1
plotFDR(bfdr1, Q=Q)

# Posterior distribution of the intercept
hist(bfdr1$coefficients[,1], 20)
abline(v=betatrue[1], col='red', lwd=4)

# Compare actual versus estimated prior probabilities of being a signal
plot(wsuccess, bfdr1$priorprob, log='xy')
abline(0,1)

# Covariate effects
plot(X[,1], log(bfdr1$priorprob/{1-bfdr1$priorprob}), ylab='Logit of prior probability')
plot(X[,2], log(bfdr1$priorprob/{1-bfdr1$priorprob}), ylab='Logit of prior probability')

# Local FDR
plot(z, bfdr1$localfdr, ylab='Local false-discovery rate')

# Extract findings at level FDR = Q
myfindings = which(bfdr1$FDR <= Q)
hist(z[myfindings], breaks=50, col='lightblue', border='blue')
table(truth = gammatrue, guess = {bfdr1$FDR <= Q})
}

}



\references{
James G. Scott, Ryan C. Kelly, Matthew A. Smith, Pengcheng Zhou, and Robert E. Kass (2015).  False discovery rate regression: application to neural synchrony detection in primary visual cortex.    Journal of the American Statistical Association, DOI: 10.1080/01621459.2014.990973. arXiv:1307.3495 [stat.ME].

Efron (2004). Large-scale simultaneous hypothesis testing: the choice of a null hypothesis. J. Amer. Statist. Assoc. 99, 96-104.

Efron (2005). Local false discovery rates. Preprint, Dept. of Statistics, Stanford University.

N.G. Polson, J.G. Scott, and J. Windle (2013.  Bayesian inference for logistic models using Polya-Gamma latent variables. Journal of the American Statistical Association (Theory and Methods) 108(504): 1339-49 (2013). arXiv:1205.0310 [stat.ME].

}

\keyword{logistic regression}
\keyword{false discovery rate}
\keyword{FDR}

% Citation:
% Adapted from <http://cran.r-project.org/doc/manuals/R-exts.html>.

% Notes:
% To check this documentation use R CMD Rdconv -t txt <filename>.

