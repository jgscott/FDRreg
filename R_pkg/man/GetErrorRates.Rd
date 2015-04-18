% package: FDRreg
\name{GetErrorRates}

\alias{GetErrorRates}

\title{Extract error rates.}

\description{
  Extract testing error rates from a fitted two-groups model.  Useful only when the ground truth is known.
}

\usage{
GetErrorRates(truth, guess)
}

\arguments{
  \item{truth}{An N dimensional binary vector indicating which entries actually are signals.}
  \item{guess}{An N dimensional binary vector indicating which entries are estimated to be signals.}
}


\value{
  \item{tpr}{The true positive rate for the guess.}
  \item{fdr}{The false discovery rate for the guess.}
  \item{confusion}{The confusion matrix.}
}

\keyword{false discovery rate}
\keyword{FDR}

% Citation:
% Adapted from <http://cran.r-project.org/doc/manuals/R-exts.html>.

% Notes:
% To check this documentation use R CMD Rdconv -t txt <filename>.

