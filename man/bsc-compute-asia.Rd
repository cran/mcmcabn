\name{bsc.compute.asia}
\docType{data}
\alias{bsc.compute.asia}

\title{Cache of pre-computed scores related to the asia dataset}
\description{This dataframe contains a cache of pre-computed scores with a maximum of two parents per node for the asia dataset.
}

     \usage{bsc.compute.asia}

\format{
     The data contains a cache of pre-computed scores with a maximum of two parents per node using \link[abn:build_score_cache]{buildscorecache} from the \code{abn} R package.
     \itemize{
     \item \code{bsc.compute.asia}: cache of score with a maximum of two parents per node.
     }}

\examples{
## This data set was generated using the following code:
library(bnlearn) #for the dataset
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

bsc.compute.asia <- buildscorecache(data.df = asia,
                                 data.dists = dist.asia,
                                 max.parents = 2)


}

\keyword{datasets}
