\name{mcmc.out.asia}
\docType{data}
\alias{mcmc.out.asia}

\title{MCMC searches from the synthetic asia dataset for use with mcmcabn library examples}
\description{This dataframe contains a cache of pre-computed scores with a maximum of two parents per node for the asia dataset.
}

     \usage{mcmc.out.asia}

\format{
     The data contains an object of class \code{mcmcabn}.
     \itemize{
     \item \code{mcmc.out.asia}: an object of class \code{mcmcabn}.
     }}

\examples{
\dontrun{
## This data set was generated using the following code:
data(asia, package='bnlearn') #for the dataset
library(abn) #for the cache of scores computing function

mcmc.out.asia <- mcmcabn(score.cache = bsc.compute.asia,
                  score = "mlik",
                  data.dists = dist.asia,
                  max.parents = 2,
                  mcmc.scheme = c(1000,99,1000),
                  seed = 42,
                  verbose = FALSE,
                  start.dag = "random",
                  prob.rev = 0.03,
                  prob.mbr = 0.03,
                  prior.choice = 2)
}
}

\keyword{datasets}
