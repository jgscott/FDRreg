% package: FDRreg
\name{FDRreg}

\alias{FDRreg}

\title{False Discovery Rate Regression}

\description{
  Estimate an empirical-Bayes false-discovery rate regression model for test statistics z and regressors X.
}


\usage{
FDRreg(z, features, nulltype='theoretical', method='pr', stderr = NULL, control=list())
}

\arguments{
  \item{z}{An N dimensional vector; z_i is the test statistic for observation i.}
  \item{features}{An N x P dimensional design matrix; feature vector x_i for the ith observation is the ith row.  This is assumed not to have a column of ones representing an intercept.  Just like in R's own lm() and glm(), the intercept will be added by the fitting algorithm.  By default the features will be centered and scaled.  This can be modified using the 'control' flag.}
  \item{nulltype}{What is the null hypothesis?  Right now the only options are Gaussian.  Choose 'empirical' for an empirical Gaussian null using Efron's central-matching method.  Choose 'theoretical' for a standard normal null (default).  Choose 'heteroscedastic' if each test statistic has a known standard error.  This last option can only be used if the 'stderr' flag is also specified.}
  \item{method}{How is the distribution of alternative test statistics estimated?  Choose 'pr' for predictive recursion (default and strongly recommended.)  For legacy reasons, there is also an 'efron' option, which uses a Poisson GLM to estimate the overall marginal density.  This is not recommended, for reasons described in the FDR regression paper.}
  \item{stderr}{Vector of standard errors.  Must be the same length as z, the vector of test statistics.  Defaults to 'NULL' in the case of a homoscedastic error model.}
  \item{control}{A list of further options to the fitting algorithm; see Details below.}
}


\value{
  \item{z}{The test statistics provided as the argument z.}
  \item{X}{The design matrix used in the regression.  This will include an added column for an intercept and will reflect the requested centering and scaling of the original features matrix.}
  \item{localfdr}{The vector of local false discovery rates (lfdr) corresponding to each element of z.  localfdr[i] is (1-p[i]). where p[i] is the fitted posterior probability that z[i] comes from the non-null (signal) population. Note localfdr is not necessarily monotonic in z, because the regression model allows the prior probability that z[i] is a signal to change with covariates X[i,].}
  \item{FDR}{The corresponding vector of cut-level false discovery rates (FDR) for the elements of z. Used for extracting findings at a given FDR level.  FDR[i] is the estimated false discovery rate for the cohort of test statistics whose local fdr's are at least as small as localfdr[i] --- that is, the z[j]'s such that localfdr[j] <= localfdr[i].  If you want the findings that meet a certain FDR cutoff, use this return value.}
  \item{x_grid}{equally-spaced midpoints of the bins used to estimate f_1(z), the density of test statistics under the alternative.}
  \item{M0}{The estimated (or assumed) null density at each of the observed z scores; M0[i] corresponds to f^0(z[i]).}
  \item{M1}{The estimated alternative density at each of the observed z scores; M1[i] corresponds to f^1(z[i]).}
  \item{fmix_grid}{The estimated value of the overall marginal density f(z) at the grid points in x_grid.}
  \item{f0_grid}{The estimated value of the null density f^0(z) at the grid points in x_grid.}
  \item{f1_grid}{The estimated value of the alternative density f^1(z) at the grid points in x_grid.}
  \item{mu0}{The mean of the null density used to fit the model.}
  \item{sig0}{The standard deviation of the null density used to fit the model.}
  \item{p0}{The estimated global null probability estimated by predictive recursion in the two-groups.  This is useful for comparison with the FDR regression results, because it assumes that the prior probability of being a signal does not change with covariates.}
 \item{priorprob}{The estimated prior probability of being a signal for each observation z[i].  Here priorprob[i] = P(z[i] is non-null). }
 \item{postprob}{The estimated posterior probabilities of being a signal each observation z[i]: postprob[i] = P(z[i] is non-null | data), and localfdr[i] = 1-postprob[i]. }
 \item{model}{The fitted FDR regression model.  A list with the coefficients, the minimized value of the complete-data negative log likelihood, and the Hessian matrix of the complete-data negative log likelihood at the mode.}

}

\details{

  This model assumes that a z-statistic z arises from

  \deqn{ f(z_i) = w_i f^1(z) + (1-w_i) f^0(z) , }

where f^1(z) and f^0(z) are the densities/marginal likelihoods under the alternative and null hypotheses, respectively, and where w_i is the prior probability that z_i is a signal (non-null case).  Predictive recursion is used to estimate f^1(z) nonparametrically; f^0(z) may either be the theoretical (standard normal) null, or an empirical null which can be estimated using the middle 25 percent of the data.  The prior probabilities w_i are estimated via logistic regression against covariates using the expectation-maximization (EM) algorithm outlined in the FDR regression companion paper.


  The control parameter can be a list with any of the following elements:
  \itemize{
    \item center, scale: Should the feature matrix be centered and scaled?  Both default to true.
    \item gridsize:The number of grid points used to estimate the alternative density.  Default: 300.
    \item decay: The learning rate in the prediction recursion algorithm for estmating f^1(z).  Analogous to the learning rate in stochastic gradient descent.  Default: -0.67.
    \item npasses: The number of passes through the data forfitting f^1(z) by prediction recursion.  Each pass sweeps through the whole z vector in a random order.  Default: 10.
    \item lambda: A small ridge-regression penalty is added to the regression model to stabilize the estimate.  Default: 0.01.  The default is a small amount of regularization when the covariates are standardized, equivalent to a Gaussian prior with variance 100.
    \item nmids, densknots: Only used if method='efron'; this functionality is retained for legacy reasons and is not recommended for actual data analysis.  nmids: How many histogram midpoints are used in Efron's spline-based estimator of f(z).  Defaults to 150.  densknots: how many spline knots.  Defaults to 10.
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

\keyword{logistic regression}
\keyword{false discovery rate}
\keyword{FDR}

% Citation:
% Adapted from <http://cran.r-project.org/doc/manuals/R-exts.html>.

% Notes:
% To check this documentation use R CMD Rdconv -t txt <filename>.

