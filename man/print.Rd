% print-mcmcabn.Rd ---
% Author           : Gilles Kratzer
% Created on :       28.02.2019
% Last modification :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{print.mcmcabn}
\alias{print.mcmcabn}
\title{Methods for mcmcabn objects}

\description{

  Method for computing on \code{mcmcabn} objects.

}

\usage{
\method{print}{mcmcabn}(x, \dots)
}

\arguments{
  \item{x}{ an object of class \code{mcmcabn}.}
  \item{\dots}{additional arguments passed to \code{print}.}
  }

\details{
There exists a \code{\link{summary}} S3 function that displays more details.
}

\author{Gilles Kratzer}

\examples{
## Example from the asia dataset from Lauritzen and Spiegelhalter (1988) provided by Scutari (2010)

print(mcmc.out.asia)
}
