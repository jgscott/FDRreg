% package: FDRreg
\name{plotFDR}

\alias{plotFDR}

\title{Plot an FDR regression model.}

\description{
  Plots the results of a fitted FDR regression model from FDRreg.
}

\details{
It is important to remember that localfdr (and therefore global FDR) is not necessarily monotonic in z, because the regression model allows the prior probability that z[i] is a signal to change with covariates x[i].
}

\usage{
plotFDR(fdrr, Q=0.1, showrug=TRUE, showfz=TRUE, showsub=TRUE) 
}

\arguments{
  \item{fdrr}{A fitted model object from FDRreg.}
  \item{Q}{The desired level at which FDR should be controlled. Defaults to 0.1, or 10 percent.}
  \item{showrug}{Logical flag indicating whether the findings at the specified FDR level should be displayed in a rug plot beneath the histogram.  Defaults to TRUE.}
  \item{showfz}{Logical flag indicating the fitted marginal density f(z) should be plotted.  Defaults to TRUE.}
  \item{showsub}{Logical flag indicating whether a subtitle should be displayed describing features of the plot.  Defaults to TRUE.}
}


\value{
No return value.
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
fdr1 <- FDRreg(z, covars=X, nmc=2500, nburn=100, nmids=120, nulltype='theoretical')
# Show the empirical-Bayes estimate of the mixture density
# and the findings at a specific FDR level
Q = 0.1
plotFDR(fdr1, Q=Q, showfz=TRUE)
}
}

%\references{
%J.G. Scott, R. Kelly, M.A. Smith, P. Zhou, and R.E. Kass (2013). False discovery rate regression: application to neural synchrony detection in primary visual cortex.  arXiv:1307.3495 [stat.ME].
%
%N.G. Polson, J.G. Scott, and J. Windle (2013.  Bayesian inference for logistic models using Polya-Gamma latent variables. Journal of the American Statistical Association (Theory and Methods) 108(504): 1339-49 (2013). arXiv:1205.0310 [stat.ME].
%
%Efron (2004). Large-scale simultaneous hypothesis testing: the choice of a null hypothesis. J. Amer. Statist. Assoc. 99, 96-104.

%Efron (2005). Local false discovery rates. Preprint, Dept. of Statistics, Stanford University.
%
%}
%
%\keyword{logistic regression}
%\keyword{false discovery rate}
%\keyword{FDR}
%
% Citation:
% Adapted from <http://cran.r-project.org/doc/manuals/R-exts.html>.
%
% Notes:
% To check this documentation use R CMD Rdconv -t txt <filename>.

