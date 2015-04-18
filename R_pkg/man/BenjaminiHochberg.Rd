% package: FDRreg
\name{BenjaminiHochberg}

\alias{BenjaminiHochberg}

\title{Run the Benjamini-Hochberg procedure.}

\description{
  Runs the Benjamini-Hochberg procedure on a vector of z statistics.
}

\usage{
BenjaminiHochberg(zscores, fdr_level)
}

\arguments{
  \item{zscores}{An N dimensional vector of z scores.}
  \item{fdr_level}{The desired level of control over the false discovery rate.}
}


\value{
  A binary vector of length N, where 0=nofinding, 1=finding at given FDR level
}

\references{
Y. Benjamini and Y. Hochberg. Controlling the false-discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society, Series B, 57: 289-300, 1995.
}


\keyword{false discovery rate}
\keyword{FDR}

% Citation:
% Adapted from <http://cran.r-project.org/doc/manuals/R-exts.html>.

% Notes:
% To check this documentation use R CMD Rdconv -t txt <filename>.

