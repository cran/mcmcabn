## ----setup, include = FALSE, cache = FALSE------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, 
  comment = "#>",
  fig.width = 7, 
  fig.height = 5,  
  fig.dpi = 200, 
  fig.align = "center",
  warning = FALSE,
  message = FALSE
)
options(digits = 3)

## ----eval=FALSE---------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install(c("RBGL", "Rgraphviz", "graph"),  version = "3.8")
#  
#  install.packages("mcmcabn")

## ----eval=TRUE----------------------------------------------------------------
library(mcmcabn)

## ---- warning = FALSE, message = FALSE----------------------------------------
library(abn) # to pre-compute the scores 
library(ggplot2) # plotting
library(ggpubr) # plotting
library(cowplot) # plotting
library(ggdag) # to plot the DAG

#plot the BN as described in the paper
dag  <- dagify(Tuberculosis~Asia,
               Either~Tuberculosis,
               XRay~Either,
               Dyspnea~Either,
               Bronchitis~Smoking,
               LungCancer~Smoking,
               Either~LungCancer,
               Dyspnea~Bronchitis)

ggdag_classic(dag, size = 6) + theme_dag_blank()


## ---- eval=FALSE--------------------------------------------------------------
#  #loglikelihood scores
#  abnCache.2par.asia <- buildScoreCache(data.df = asia,
#                                      data.dists = dist.asia,
#                                      max.parents = 2)

## ---- eval=FALSE--------------------------------------------------------------
#  mcmc.2par.asia <- mcmcabn(score.cache = abnCache.2par.asia,
#                           score = "mlik",
#                           data.dists = dist.asia,
#                           max.parents = 2,
#                           mcmc.scheme = c(1000,99,1000),
#                           seed = 42,
#                           verbose = FALSE,
#                           start.dag = "random",
#                           prob.rev = 0.03,
#                           prob.mbr = 0.03,
#                           prior.choice = 2)

## ---- echo=FALSE--------------------------------------------------------------
##to speed up building of the vignette, we store the output
data("mcmc_run_asia")

## -----------------------------------------------------------------------------
#maximum score get from the MCMC sampler
max(mcmc.2par.asia$scores)

#maximum scoring network using exact search (not MCMC based) 
dag <- mostProbable(score.cache = abnCache.2par.asia)
fitAbn(object = dag)$mlik

## -----------------------------------------------------------------------------
plot(mcmc.2par.asia)

## -----------------------------------------------------------------------------
plot(mcmc.2par.asia, max.score = TRUE)

## -----------------------------------------------------------------------------
summary(mcmc.2par.asia)

## -----------------------------------------------------------------------------
query(mcmcabn = mcmc.2par.asia)

## -----------------------------------------------------------------------------
query(mcmcabn = mcmc.2par.asia, formula = ~ LungCancer|Smoking)

## -----------------------------------------------------------------------------
query(mcmcabn = mcmc.2par.asia, formula = ~ LungCancer|Smoking + Bronchitis|Smoking)

## -----------------------------------------------------------------------------
query(mcmcabn = mcmc.2par.asia ,formula = ~LungCancer|Smoking + Bronchitis|Smoking - Tuberculosis|Smoking - XRay|Bronchitis)

