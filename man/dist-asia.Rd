\name{asia}
\docType{data}
\alias{dist.asia}
\alias{asia}

\title{The \code{asia} dataset}
\description{The \code{asia} dataset and named list of distribution to analyze asia dataset.
}

\usage{
asia
dist.asia}
\details{
The data is a synthetic dataset from Lauritzen and Spiegelhalter (1988) provided by Scutari (2010).
Renamed for convenience, similar to the use in the package \pkg{abn}.
}
\format{
     \itemize{
     \item \code{asia}: a dataframe with 500 observations and 8 variables.
     Each variable is encoded as a factor with levels \code{"no"} and \code{"yes"}.
     \item \code{dist.asia}: a named list giving the distribution for each node in the network (here \code{"binomial"} for all).
}}

\examples{
\dontrun{
## This data set was generated using the following code:
data(asia, package = 'bnlearn') # for the dataset

# Renaming columns of the dataset
colnames(asia) <- c("Asia", "Smoking", "Tuberculosis", "LungCancer",
  "Bronchitis", "Either", "XRay", "Dyspnea")

# Defining the distribution list
dist.asia <- list(Asia = "binomial", Smoking = "binomial",
    Tuberculosis = "binomial", LungCancer = "binomial",
    Bronchitis = "binomial", Either = "binomial",
    XRay = "binomial", Dyspnea = "binomial")
}
}
\keyword{datasets}
