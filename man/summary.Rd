% summary.mcmcabn.Rd ---
% Author           : Gilles Kratzer
% Created on :       28.02.2019
% Last modification :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{summary.mcmcabn}
\alias{summary.mcmcabn}
\title{Function to summarize MCMC run generated by mcmcabn}

\usage{
\method{summary}{mcmcabn}(object,
     quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
     lag.max = 10,
     \dots)
     }

\arguments{
  \item{object}{object of class \code{mcmcabn}.}
  \item{quantiles}{numeric vector of probabilities with values in [0,1]. (Values up to 2e-14 outside that range are accepted and moved to the nearby endpoint.)}
  \item{lag.max}{maximum lag at which to calculate the \link{acf}. Default is set to 10.}
  \item{\dots}{arguments to be passed to methods.}
  }

\description{Summary method for mcmcabn objects.}

\details{The summary function for \code{mcmcabn} objects returns multiple summary metrics for assessing the quality of the MCMC run. Thinning is the number of thinned MCMC steps for one MCMC returned.}

\value{This method prints: the number of burn-in steps, the number of MCMC steps, the thinning, the maximum achieved score, the empirical mean of the MCMC samples, the empirical standard deviation of the MCMC samples, the user defined quantiles of the posterior network score, the global acceptance rate, a table of the accepted and rejected moves in function of the methods used, the sample size adjusted for autocorrelation and the autocorrelations by lag.}

\author{Gilles Kratzer}

\references{
Scutari, M. (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1 - 22. doi:http://dx.doi.org/10.18637/jss.v035.i03.
}


\examples{
## Example from the asia dataset from Lauritzen and Spiegelhalter (1988)
## provided by Scutari (2010)
#summary the MCMC run
summary(mcmc.2par.asia)
}
