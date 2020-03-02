% print-summar.mcmcabn.Rd ---
% Author           : Gilles Kratzer
% Created on :       31.10.2019
% Last modification :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{print.summary.mcmcabn}
\alias{print.summary.mcmcabn}
\title{Methods for printing the summary of mcmcabn objects}

\description{

  Method for printing the summary of \code{mcmcabn} objects.

}

\usage{
  \method{print}{summary.mcmcabn}(x, \dots)
}

\arguments{
  \item{x}{an object of class \code{summary.mcmcabn}.}
  \item{\dots}{additional arguments passed to \code{print}.}
  }

\details{
There exists a \code{\link{summary}} S3 function that displays more details.
}

\author{Gilles Kratzer}

\examples{
## Example from the asia dataset from Lauritzen and Spiegelhalter (1988)
## provided by Scutari (2010)

summary(mcmc.out.asia)
}
