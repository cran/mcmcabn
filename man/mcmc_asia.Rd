\name{mcmc_asia}
\docType{data}
\alias{mcmc_run_asia}
\alias{mcmc.run.asia}
\alias{mcmc.4par.asia}
\alias{mcmc.2par.asia}

\title{Precomputed \code{mcmcabn} objects.}

\description{10^5 MCMC runs with 1000 burn-in runs from the synthetic
  \code{asia} synthetic dataset. Used in examples and vignettes of the package \pkg{mcmcabn}.
}

\usage{
mcmc.2par.asia
mcmc.4par.asia
}

\format{
     The objects are of class \code{mcmcabn} and rely on a cache of score computed using
     \link[abn:buildScoreCache]{buildScoreCache} from the R package \pkg{abn}.
}

\seealso{\code{\link{abnCache.4par.asia}},\code{\link{asia}}}

\examples{
## Compare the scores based on 4 and 2 parents:
plot(density(mcmc.4par.asia$scores))
lines(density(mcmc.2par.asia$scores), col = 2)

## Best scores do not match dags with 4 parents:
plot(mcmc.4par.asia$scores, type='l', col='gray')
maxpar <- apply(mcmc.4par.asia$dags, 3, function(x) max(rowSums(x)))
rug((1:1001)[maxpar==4], side = 3, col = 2)


\dontrun{
## This data set was generated using the following code:
## (The following lines take some time to evaluate).
require(abn)
set.seed(3)
abnCache.4par.asia <- buildScoreCache(data.df = asia,
    data.dists = dist.asia, max.parents = 4)

mcmc.4par.asia <- mcmcabn(score.cache = abnCache.4par.asia,
    score = "mlik", data.dists = dist.asia, max.parents = 4,
    mcmc.scheme = c(1000, 99, 1000), seed = 42, verbose = FALSE,
    start.dag = "random", prob.rev = 0.03, prob.mbr = 0.03,
    prior.choice = 2)


set.seed(3)
abnCache.2par.asia <- buildScoreCache(data.df = asia,
    data.dists = dist.asia, max.parents = 2)

mcmc.2par.asia <- mcmcabn(score.cache = abnCache.2par.asia,
    score = "mlik", data.dists = dist.asia, max.parents = 2,
    mcmc.scheme = c(1000, 99, 1000), seed = 42, verbose = FALSE,
    start.dag = "random", prob.rev = 0.03, prob.mbr = 0.03,
    prior.choice = 2)

}


}

\keyword{datasets}
