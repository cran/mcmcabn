\name{abnCache.asia}
\docType{data}
\alias{bsc.compute.asia}
\alias{abnCache.asia}
\alias{abnCache.2par.asia}
\alias{abnCache.4par.asia}

\title{Cache of pre-computed scores related to the asia dataset}
\description{This list contains a cache of pre-computed scores with a maximum of two/four parents per node for the asia dataset.
}

\usage{
abnCache.2par.asia
abnCache.4par.asia
}

\format{
     The lists are objects of class \code{abnCache} and contain a cache of pre-computed scores with a maximum of two parents (of four parents) per node using \link[abn:build_score_cache]{buildScoreCache} from the \code{abn} R package.
     \itemize{
     \item \code{abnCache.2par.asia}: list of 14 containing a cache of score with a maximum of two parents per node.
     \item \code{abnCache.4par.asia}: list of 14 containing a cache of score with a maximum of four parents per node.
 }}

\seealso{\code{\link{mcmc.2par.asia}},\code{\link{asia}}}
\examples{
## The object was generated using the following code:
\dontrun{
require("abn")
abnCache.2par.asia <- buildScoreCache(data.df = asia,
    data.dists = dist.asia, max.parents = 2)

abnCache.4par.asia <- buildScoreCache(data.df = asia,
    data.dists = dist.asia, max.parents = 4)
}}

\keyword{datasets}
