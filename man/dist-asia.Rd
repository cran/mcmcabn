\name{dist.asia}
\docType{data}
\alias{dist.asia}

\title{Named list of distributions to analyze asia dataset}
\description{Named list of distribution to analyze asia dataset.
}

     \usage{dist.asia}

\format{
     The data contains a cache of pre-computed scores with a maximum of two parents per node.
     \itemize{
     \item \code{dist.asia}: a named list giving the distribution for each node in the network.
     }}

\examples{
## This data set was generated using the following code:
library(bnlearn) #for the dataset

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
}

\keyword{datasets}
