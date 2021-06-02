\name{mcmc_run_asia}
\docType{data}
\alias{mcmc_run_asia}

\title{List of files to reproduce examples \code{mcmcabn} library.}

\description{10^5 MCMC runs with 1000 burn-in runs from the asia synthetic dataset from Lauritzen and Spiegelhalter (1988) provided by Scutari (2010). A named list of distributions and pre-computed scores.
}

     \usage{data("mcmc_run_asia")}

\format{
     The data contains an object of class mcmcabn and a cache of score computed using
     \link[abn:build_score_cache]{buildScoreCache} from the \code{abn} R package.
     \itemize{
     \item \code{bsc.compute.asia}: cache of score with a maximum of two parents per node computed using \link[abn:build_score_cache]{buildScoreCache} from the \code{abn} R package.
     \item \code{dist.asia}: a named list giving the distribution for each node in the network.
     \item \code{mcmc.out.asia}: an object of class mcmcabn.
     }}

\examples{
\dontrun{
## This data set was generated using the following code:
data(asia, package='bnlearn') #for the dataset
library(abn) #for the cache of score function

#renaming columns of the dataset
colnames(asia) <- c("Asia",
  "Smoking",
  "Tuberculosis",
  "LungCancer",
  "Bronchitis",
  "Either",
  "XRay",
  "Dyspnea")

#lets define the distribution list
dist.asia <- list(Asia = "binomial",
            Smoking = "binomial",
            Tuberculosis = "binomial",
            LungCancer = "binomial",
            Bronchitis = "binomial",
            Either = "binomial",
            XRay = "binomial",
            Dyspnea = "binomial")

bsc.compute.asia <- buildScoreCache(data.df = asia,
                                 data.dists = dist.asia,
                                 max.parents = 2)

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
