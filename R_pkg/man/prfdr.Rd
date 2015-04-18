% package: FDRreg
\name{prfdr}

\alias{prfdr}

\title{Predictive recursion for the two-groups model.}

\description{
  Estimate the two-groups model by predictive recursion.
}

\usage{
prfdr(z, mu0=0.0, sig0=1.0, control=list())
}

\arguments{
  \item{z}{An N dimensional vector; z_i is the test statistic for observation i.}
   \item{mu0}{The mean of the null density used to fit the model.}
  \item{sig0}{The standard deviation of the null density used to fit the model.}
  \item{control}{A list of further options to the fitting algorithm; see Details below.}
}


\value{
  \item{x_grid}{Equally-spaced midpoints of the bins used to estimate f_1(z), the density of test statistics under the alternative.}
  \item{pi0}{The estimated global null probability (1-w).}
  \item{pitheta_grid}{The estimated value of the prior pi(theta)}
  \item{f1_grid}{The estimated value of the alternative density f^1(z) at the grid points in x_grid.}
  \item{f0_grid}{The estimated value of the null density f^0(z) at the grid points in x_grid.}
  \item{fmix_grid}{The estimated value of the overall marginal density f(z) at the grid points in x_grid.}


  \item{f0_z}{The null density at each of the observed z scores.}
  \item{f1_z}{The estimated alternative density at each of the observed z scores.}
  \item{fmix_z}{The estimated overall mixture density at each of the observed z scores.}

  \item{mu0}{The mean of the null density used to fit the model.}
  \item{sig0}{The standard deviation of the null density used to fit the model.}
 \item{postprob}{The estimated posterior probabilities of being a signal each observation z[i]: postprob[i] = P(z[i] is non-null | data), and localfdr[i] = 1-postprob[i]. }

}

\details{


  This model assumes that a z-statistic z arises from

  \deqn{ f(z_i) = w f^1(z) + (1-w) f^0(z) , }

where f^1(z) and f^0(z) are the densities/marginal likelihoods under the alternative and null hypotheses, respectively, and where w is the prior probability that z_i is a signal (non-null case).  Predictive recursion is used to estimate f^1(z) nonparametrically; f^0(z) may either be the theoretical (standard normal) null, or an empirical null which can be estimated using the middle 25 percent of the data. 

The assumption of predictive recursion is that f^1(z) is a convolution of the Gaussian null density with respect to a prior distribution pi(theta), which is estimated by the method.


  The control parameter can be a list with any of the following elements:
  \itemize{
    \item gridsize: The number of grid points used to estimate the alternative density.  Default: 500.
    \item decay: The learning rate in the prediction recursion algorithm for estimating f^1(z).  Analogous to the learning rate in stochastic gradient descent.  Default: -0.67.
    \item npasses: The number of passes through the data for fitting f^1(z) by prediction recursion.  Each pass sweeps through the whole z vector in a random order.  Default: 10.
  }
}

\examples{

library(FDRreg)

# Simulated data
N = 10000
wsuccess = 0.05

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

# Fit the model by PR and show the fit
pr1 = prfdr(z)
hist(z, 100, col='lightblue', border=NA, prob=TRUE)
lines(pr1$x_grid, pr1$fmix_grid, col='blue')

# Show estimated prior versus actual prior
plot(pr1$x_grid, pr1$pitheta_grid, col='red', type='l')
curve(dnorm(x, 3, 0.5), add=TRUE)

}
}

\references{
James G. Scott, Ryan C. Kelly, Matthew A. Smith, Pengcheng Zhou, and Robert E. Kass (2015).  False discovery rate regression: application to neural synchrony detection in primary visual cortex.    Journal of the American Statistical Association, DOI: 10.1080/01621459.2014.990973. arXiv:1307.3495 [stat.ME].

R. Martin and S. Tokdar. A nonparametric empirical Bayes framework for large-scale multiple testing. Biostatistics, 13(3):427-39, 2012.

}

\keyword{predictive recursion}
\keyword{nonparametric empirical Bayes}
\keyword{false discovery rate}
\keyword{FDR}

% Citation:
% Adapted from <http://cran.r-project.org/doc/manuals/R-exts.html>.

% Notes:
% To check this documentation use R CMD Rdconv -t txt <filename>.

