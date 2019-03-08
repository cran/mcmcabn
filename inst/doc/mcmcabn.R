## ----setup, include = FALSE, cache = FALSE-------------------------------
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

## ----eval=FALSE----------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install(c("RBGL", "Rgraphviz", "graph"),  version = "3.8")
#  
#  install.packages("mcmcabn")

## ----eval=TRUE-----------------------------------------------------------
library(mcmcabn)

## ---- warning = FALSE, message = FALSE-----------------------------------
library(bnlearn) #for the dataset
library(abn) #to pre-compute the scores 
library(ggplot2) #plotting
library(ggpubr) #plotting
library(cowplot) #plotting
library(ggdag) #to plot the DAG

#renaming columns of the dataset
colnames(asia) <- c("Asia", "Smoking", "Tuberculosis", "LungCancer", "Bronchitis", "Either", "XRay", "Dyspnea")

#lets define the distribution list
dist.asia <- list(Asia = "binomial",
                  Smoking = "binomial",
                  Tuberculosis = "binomial",
                  LungCancer = "binomial", 
                  Bronchitis = "binomial",
                  Either = "binomial",
                  XRay = "binomial", 
                  Dyspnea = "binomial")

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


## ---- eval=FALSE---------------------------------------------------------
#  #loglikelihood scores
#  bsc.compute.asia <- buildscorecache(data.df = asia,
#                                      data.dists = dist.asia,
#                                      max.parents = 2)

## ---- eval=FALSE---------------------------------------------------------
#  mcmc.out.asia <- mcmcabn(score.cache = bsc.compute.asia,
#                           score = "mlik",
#                           data.dists = dist.asia,
#                           max.parents = 2,
#                           mcmc.scheme = c(1000,99,10000),
#                           seed = 42,
#                           verbose = FALSE,
#                           start.dag = "random",
#                           prob.rev = 0.03,
#                           prob.mbr = 0.03,
#                           prior.choice = 2)

## ---- echo=FALSE---------------------------------------------------------
##to speed up building of the vignette, we store the output
data("mcmc_run_asia")

## ------------------------------------------------------------------------
#maximum score get from the MCMC sampler
max(mcmc.out.asia$scores)

#maximum scoring network using exact search (not MCMC based) 
dag <- mostprobable(score.cache = bsc.compute.asia)
fitabn(dag.m = dag,data.df = asia, data.dists = dist.asia)$mlik

## ------------------------------------------------------------------------
plot(mcmc.out.asia)

## ------------------------------------------------------------------------
plot(mcmc.out.asia, max.score = TRUE)

## ------------------------------------------------------------------------
summary(mcmc.out.asia)

## ------------------------------------------------------------------------
query(mcmcabn = mcmc.out.asia)

## ------------------------------------------------------------------------
query(mcmcabn = mcmc.out.asia, formula = ~ LungCancer|Smoking)

## ------------------------------------------------------------------------
query(mcmcabn = mcmc.out.asia, formula = ~ LungCancer|Smoking + Bronchitis|Smoking)

## ------------------------------------------------------------------------
query(mcmcabn = mcmc.out.asia ,formula = ~LungCancer|Smoking + Bronchitis|Smoking - Tuberculosis|Smoking - XRay|Bronchitis)

