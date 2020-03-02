
\name{. mcmcabn .}
\alias{overview}
\alias{mcmcabn-package}
\docType{package}
\title{mcmcabn Package}
\description{\code{mcmcabn} is a structural MCMC sampler for Directed Acyclic Graphs (DAGs). The three main problems that can be addressed by this R package are selecting the most probable structure based on a cache of pre-computed scores, controlling for overfitting and sampling the landscape of high scoring structures. It allows to quantify the marginal impact of relationships of interest by marginalising out over structures or nuisance dependencies.
}

\section{General overview}{What is \pkg{mcmcabn}:\cr

Bayesian network modeling is a data analysis technique that is ideally suited to messy, highly correlated, and complex datasets. This methodology is somewhat distinct from other forms of statistical modeling in that its focus is on structure discovery - determining an optimal graphical model that describes the inter-relationships in the underlying processes which generated the data. It is a multivariate technique and can used for one or many dependent variables. This is a data-driven approach, as opposed to, rely only on subjective expert opinion to determine how variables of interest are inter-related (for example, structural equation modeling).

The R package mcmcabn is a structural MCMC sampler for Directed Acyclic Graphs (DAGs). It contains routines to compute, analyze, and report MCMC samples. This structural sampler supports the new edge reversal move from Grzegorczyk and Husmeier (2008) <doi:10.1007/s10994-008-5057-7> and the Markov blanket resampling from Su and Borsuk (2016) <http://jmlr.org/papers/v17/su16a.html>. It supports three priors: a prior controlling for structure complexity from Koivisto and Sood (2004) <http://dl.acm.org/citation.cfm?id=1005332.1005352>, an uninformative prior and a user-defined prior.

The web pages \url{https://www.math.uzh.ch/pages/mcmcabn/} provide further case studies. See also the files stored in the package directories
\code{inst/FCV_code}.

}
\references{


Kratzer G, Lewis FI, Willi B, Meli ML, Boretti FS, Hofmann-Lehmann R, Torgerson P, Furrer R and Hartnack S (2020) Bayesian Network Modeling Applied to Feline Calicivirus Infection Among Cats in Switzerland. Front. Vet. Sci. 7:73. doi: 10.3389/fvets.2020.00073 (Preferred reference)

Kratzer, G. and Furrer, R. (2019). Is a single unique Bayesian network enough to accurately represent
  your data?. arXiv preprint arXiv:1902.06641

  Kratzer, G. and Furrer, R. (2019). mcmcabn: a structural MCMC sampler for DAGs learned from observed
  systemic datasets. R package version 0.3.  https://CRAN.R-project.org/package=mcmcabn
  }


\examples{
## Citations:
citation('mcmcabn')
}
\author{Gilles Kratzer}
\keyword{documentation}
\keyword{package}
