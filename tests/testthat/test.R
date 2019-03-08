###############################################################################
## test.R ---
## Author          : Gilles Kratzer
## Document created: 25/02/2019
###############################################################################

##Purpose: Test the mcmcabn software

Sys.setenv("R_TESTS" = "")
#library(testthat)
#library(abn)
library(bnlearn)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##General tests
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

context("General tests")

test_that("mc3",{

  #add arrow
  dist <- list(A="binomial",
               B="binomial",
               C="binomial",
               D="binomial",
               E="binomial")

  data.param <- matrix(data = c(0,0,0.9,0,0,0.9,0,0.9,0,0.9,0,0,0,0,0,0,0,0,0,0.9,0,0,0.9,0,0.9),nrow = 5L,ncol = 5L,byrow = T)

  diag(data.param)<-0.5
  colnames(data.param) <- rownames(data.param) <- names(dist)

  out.sim.0 <- simulateabn(data.dists = dist,n.chains = 1,n.adapt = 1000,n.thin = 1,n.iter = 1000,data.param = 0.4*data.param, simulate = TRUE,seed = 132,verbose = FALSE)

  bsc.compute.0 <- buildscorecache(data.df = out.sim.0, data.dists = dist,  max.parents = 3)

  mc3.out <- mcmcabn(score.cache = bsc.compute.0,
                  score = "mlik",
                  data.dists = dist,
                  max.parents = 3,
                  mcmc.scheme = c(1000,0,0),
                  seed = 42,
                  verbose = FALSE,
                  start.dag = "random",
                  prob.rev = 0,
                  prob.mbr = 0,
                  prior.choice = 1)

  dag <- mostprobable(score.cache = bsc.compute.0,prior.choice = 1)
  expect_equal(max(mc3.out$scores),fitabn(dag.m = dag,data.df = out.sim.0,data.dists = dist)$mlik)

  mc3.out <- mcmcabn(score.cache = bsc.compute.0,
                     score = "mlik",
                     data.dists = dist,
                     max.parents = 3,
                     mcmc.scheme = c(5000,0,0),
                     seed = 56,
                     verbose = FALSE,
                     start.dag = "random",
                     prob.rev = 0,
                     prob.mbr = 0,
                     prior.choice = 2)

  dag <- mostprobable(score.cache = bsc.compute.0,prior.choice = 2)
  expect_equal(round(max(mc3.out$scores),digits = 0),round(fitabn(dag.m = dag,data.df = out.sim.0,data.dists = dist)$mlik,digits = 0))

  data.param.eq <- matrix(data = 0,nrow = 5,ncol = 5)
  mc3.out <- mcmcabn(score.cache = bsc.compute.0,
                     score = "mlik",
                     data.dists = dist,
                     max.parents = 1,
                     mcmc.scheme = c(1000,0,0),
                     seed = 42,
                     verbose = FALSE,
                     start.dag = "random",
                     prob.rev = 0,
                     prob.mbr = 0,
                     prior.dag = data.param.eq,
                     prior.lambda = 5,
                     prior.choice = 3)

  expect_false(table(apply(mc3.out$dags,3, sum))[1]<500)

  data.param.eq <- matrix(data = 1,nrow = 5,ncol = 5)
  mc3.out <- mcmcabn(score.cache = bsc.compute.0,
                     score = "mlik",
                     data.dists = dist,
                     max.parents = 1,
                     mcmc.scheme = c(1000,0,0),
                     seed = 42,
                     verbose = FALSE,
                     start.dag = "random",
                     prob.rev = 0,
                     prob.mbr = 0,
                     prior.dag = data.param.eq,
                     prior.lambda = 5,
                     prior.choice = 3)

  expect_false(table(apply(mc3.out$dags,3, sum))[1]>500)

})

test_that("REV",{

  dist <- list(A="binomial",
               B="binomial",
               C="binomial",
               D="binomial",
               E="binomial")

  data.param <- matrix(data = c(0,0.9,0,0.9,0,0,0,0,0,0,0.9,0.9,0,0,0,0,0,0,0.9,0,0,0,0,0.9,0),nrow = 5L,ncol = 5L,byrow = T)

  diag(data.param)<-0.5
  colnames(data.param) <- rownames(data.param) <- names(dist)

  out.sim.0 <- invisible(simulateabn(data.dists = dist,n.chains = 1,n.adapt = 1000,n.thin = 1,n.iter = 1000,data.param = 0.4*data.param, simulate = TRUE,seed = 132))

  bsc.compute.0 <- buildscorecache(data.df = out.sim.0, data.dists = dist,  max.parents = 3)

  mc3.out <- mcmcabn(score.cache = bsc.compute.0,
                     score = "mlik",
                     data.dists = dist,
                     max.parents = 3,
                     mcmc.scheme = c(100,0,0),
                     seed = 23213,
                     verbose = FALSE,
                     start.dag = "random",
                     prob.rev = 0.5,
                     prob.mbr = 0,
                     prior.choice = 1)

  dag <- mostprobable(score.cache = bsc.compute.0,prior.choice = 1)
  expect_equal(max(mc3.out$scores),fitabn(dag.m = dag,data.df = out.sim.0,data.dists = dist)$mlik)

  expect_silent(mcmcabn(score.cache = bsc.compute.0,
                        score = "mlik",
                        data.dists = dist,
                        max.parents = 3,
                        mcmc.scheme = c(100,0,0),
                        seed = 32132,
                        verbose = FALSE,
                        start.dag = "random",
                        prob.rev = 1,
                        prob.mbr = 0,
                        prior.choice = 1))

})

test_that("MBR",{

  dist <- list(A="binomial",
               B="binomial",
               C="binomial",
               D="binomial",
               E="binomial")

  data.param <- matrix(data = c(0,0,0.9,0.9,0,0,0.9,0,0,0.9,0,0,0,0,0,0,0,0,0,0.9,0,0,0,0,0),nrow = 5L,ncol = 5L,byrow = T)

  diag(data.param)<-0.5
  colnames(data.param) <- rownames(data.param) <- names(dist)

  out.sim.0 <- invisible(simulateabn(data.dists = dist,n.chains = 1,n.adapt = 1000,n.thin = 1,n.iter = 1000,data.param = 0.4*data.param, simulate = TRUE,seed = 132))

  bsc.compute.0 <- buildscorecache(data.df = out.sim.0, data.dists = dist,  max.parents = 3)

  mc3.out <- mcmcabn(score.cache = bsc.compute.0,
                     score = "mlik",
                     data.dists = dist,
                     max.parents = 3,
                     mcmc.scheme = c(100,0,0),
                     seed = 23213,
                     verbose = FALSE,
                     start.dag = "random",
                     prob.rev = 0,
                     prob.mbr = 0.1,
                     prior.choice = 2)

  dag <- mostprobable(score.cache = bsc.compute.0,prior.choice = 1)
  expect_equal(max(mc3.out$scores),fitabn(dag.m = dag,data.df = out.sim.0,data.dists = dist)$mlik)

  expect_silent(mcmcabn(score.cache = bsc.compute.0,
                        score = "mlik",
                        data.dists = dist,
                        max.parents = 3,
                        mcmc.scheme = c(100,0,0),
                        seed = 32132,
                        verbose = FALSE,
                        start.dag = "random",
                        prob.rev = 0,
                        prob.mbr = 1,
                        prior.choice = 1))

})

test_that("mcmcabn",{

  dist<-list(A="binomial",
             B="binomial",
             C="binomial",
             D="binomial",
             E="binomial")

  data.param.0 <- matrix(data = c(0,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1),nrow = 5L,ncol = 5L,byrow = T)
  colnames(data.param.0) <- rownames(data.param.0) <- names(dist)

  out.sim.0 <- invisible(simulateabn(data.dists = dist,n.chains = 1,n.adapt = 20,n.thin = 1,n.iter = 100,data.param = data.param.0, simulate = TRUE,seed = 132))

  bsc.compute.0 <- buildscorecache(data.df = out.sim.0, data.dists = dist,  max.parents = 2)

  expect_error(mcmcabn(score.cache = bsc.compute.0,
                       score = "blabla",
                       data.dists = dist,
                       max.parents = 1,
                       mcmc.scheme = c(5,0,0),
                       seed = 2343,
                       verbose = FALSE,
                       start.dag = "random",
                       prob.rev = 0,
                       prob.mbr = 0,
                       prior.choice = 1))

  expect_error(mcmcabn(score.cache = bsc.compute.0,
                       score = "mlik",
                       data.dists = unname(dist),
                       max.parents = 1,
                       mcmc.scheme = c(5,0,0),
                       seed = 2343,
                       verbose = FALSE,
                       start.dag = "random",
                       prob.rev = 0,
                       prob.mbr = 0,
                       prior.choice = 1))

  expect_error(mcmcabn(score.cache = bsc.compute.0,
                       score = "mlik",
                       data.dists = unname(dist),
                       max.parents = 0,
                       mcmc.scheme = c(5,0,0),
                       seed = 2343,
                       verbose = FALSE,
                       start.dag = "random",
                       prob.rev = 0,
                       prob.mbr = 0,
                       prior.choice = 1))

  expect_error(mcmcabn(score.cache = bsc.compute.0,
                       score = "mlik",
                       data.dists = dist,
                       max.parents = 0,
                       mcmc.scheme = c(-1,0,0),
                       seed = 2343,
                       verbose = FALSE,
                       start.dag = "random",
                       prob.rev = 0,
                       prob.mbr = 0,
                       prior.choice = 1))

  expect_error(mcmcabn(score.cache = bsc.compute.0,
                       score = "mlik",
                       data.dists = dist,
                       max.parents = 1,
                       mcmc.scheme = c(1,0,0),
                       seed = 2343,
                       verbose = FALSE,
                       start.dag = "random",
                       prob.rev = 0,
                       prob.mbr = 0,
                       prior.choice = 15))

  expect_error(mcmcabn(score.cache = bsc.compute.0,
                       score = "mlik",
                       data.dists = dist,
                       max.parents = 1,
                       mcmc.scheme = c(1,0,0),
                       seed = 2343,
                       verbose = FALSE,
                       start.dag = "random",
                       prob.rev = -0.1,
                       prob.mbr = 0,
                       prior.choice = 1))


  expect_error(mcmcabn(score.cache = bsc.compute.0,
                       score = "mlik",
                       data.dists = dist,
                       max.parents = 1,
                       mcmc.scheme = c(1,0,0),
                       seed = 2343,
                       verbose = FALSE,
                       start.dag = "random",
                       prob.rev = -0.1,
                       prob.mbr = 0,
                       prior.choice = 1))

##asia
  data(asia)

  dist.asia <- list(Asia = "binomial",
                    Smoking = "binomial",
                    Tuberculosis = "binomial",
                    LungCancer = "binomial",
                    Bronchitis = "binomial",
                    Either = "binomial",
                    XRay = "binomial",
                    Dyspnea = "binomial")

  colnames(asia) <- c("Asia","Smoking", "Tuberculosis", "LungCancer", "Bronchitis", "Either", "XRay", "Dyspnea")


  bsc.compute.asia <- buildscorecache(data.df = asia,
                                      data.dists = dist.asia,
                                      max.parents = 2)


  mcmc.out.asia <- mcmcabn(score.cache = bsc.compute.asia,
                           score = "mlik",
                           data.dists = dist.asia,
                           max.parents = 2,
                           mcmc.scheme = c(500,0,0),
                           seed = 456,
                           verbose = FALSE,
                           start.dag = "random",
                           prob.rev = 0.1,
                           prob.mbr = 0.1,
                           prior.choice = 2)


  #maximum scoring network using exact search (not MCMC based)
  dag <- mostprobable(score.cache = bsc.compute.asia)
  expect_equal(max(mcmc.out.asia$scores),fitabn(dag.m = dag,data.df = asia,data.dists = dist.asia)$mlik)

  mcmc.out.asia <- mcmcabn(score.cache = bsc.compute.asia,
                           score = "mlik",
                           data.dists = dist.asia,
                           max.parents = 2,
                           mcmc.scheme = c(250,0,0),
                           seed = 789,
                           verbose = FALSE,
                           start.dag = "hc",
                           prob.rev = 0.1,
                           prob.mbr = 0.1,
                           prior.choice = 2)

  expect_equal(max(mcmc.out.asia$scores),fitabn(dag.m = dag,data.df = asia,data.dists = dist.asia)$mlik)

  ## marks datasets

  data(marks)

  dist.marks <- list(MECH = "gaussian",
                    VECT = "gaussian",
                    ALG = "gaussian",
                    ANL = "gaussian",
                    STAT = "gaussian")

  #colnames(asia) <- c("Asia","Smoking", "Tuberculosis", "LungCancer", "Bronchitis", "Either", "XRay", "Dyspnea")


  bsc.compute.marks <- buildscorecache(data.df = marks,
                                      data.dists = dist.marks,
                                      max.parents = 2)


  mcmc.out.marks <- mcmcabn(score.cache = bsc.compute.marks,
                           score = "mlik",
                           data.dists = dist.marks,
                           max.parents = 2,
                           mcmc.scheme = c(250,0,0),
                           seed = 789,
                           verbose = FALSE,
                           start.dag = "random",
                           prob.rev = 0.03,
                           prob.mbr = 0.03,
                           prior.choice = 2)



  #maximum scoring network using exact search (not MCMC based)
  dag <- mostprobable(score.cache = bsc.compute.marks)
  expect_equal(max(mcmc.out.marks$scores),fitabn(dag.m = dag,data.df = marks,data.dists = dist.marks)$mlik)

  ##tests
  data(gaussian.test)

  dist.gaussian.test <- list(A = "gaussian",
                     B = "gaussian",
                     C = "gaussian",
                     D = "gaussian",
                     E = "gaussian",
                     G = "gaussian",
                     H ="gaussian")

  colnames(gaussian.test) <- c("A","B","C","D","E","G","H")

  bsc.compute.gaussian.test <- buildscorecache(data.df = gaussian.test,
                                       data.dists = dist.gaussian.test,
                                       max.parents = 2)


  mcmc.out.gaussian.test <- mcmcabn(score.cache = bsc.compute.gaussian.test,
                            score = "mlik",
                            data.dists = dist.gaussian.test,
                            max.parents = 2,
                            mcmc.scheme = c(250,0,0),
                            seed = 678,
                            verbose = FALSE,
                            start.dag = "random",
                            prob.rev = 0.1,
                            prob.mbr = 0.1,
                            prior.choice = 2)



  #maximum scoring network using exact search (not MCMC based)
  dag <- mostprobable(score.cache = bsc.compute.gaussian.test)
  expect_equal(max(mcmc.out.gaussian.test$scores),fitabn(dag.m = dag,data.df = gaussian.test,data.dists = dist.gaussian.test)$mlik)

  })

test_that("query",{

  data(gaussian.test)

  dist.gaussian.test <- list(A = "gaussian",
                             B = "gaussian",
                             C = "gaussian",
                             D = "gaussian",
                             E = "gaussian",
                             G = "gaussian",
                             H ="gaussian")

  colnames(gaussian.test) <- c("A","B","C","D","E","G","H")

  bsc.compute.gaussian.test <- buildscorecache(data.df = gaussian.test,
                                               data.dists = dist.gaussian.test,
                                               max.parents = 2)


  mcmc.out.gaussian.test <- mcmcabn(score.cache = bsc.compute.gaussian.test,
                                    score = "mlik",
                                    data.dists = dist.gaussian.test,
                                    max.parents = 2,
                                    mcmc.scheme = c(100,0,0),
                                    seed = 623178,
                                    verbose = FALSE,
                                    start.dag = "random",
                                    prob.rev = 0.1,
                                    prob.mbr = 0.1,
                                    prior.choice = 2)

expect_true(is.matrix(query(mcmcabn = mcmc.out.gaussian.test)))

expect_equal(query(mcmcabn = mcmc.out.gaussian.test)[1,2],query(mcmcabn = mcmc.out.gaussian.test,formula = ~A|B))
expect_equal(query(mcmcabn = mcmc.out.gaussian.test)[2,3],query(mcmcabn = mcmc.out.gaussian.test,formula = ~B|C))
expect_equal(query(mcmcabn = mcmc.out.gaussian.test)[6,1],query(mcmcabn = mcmc.out.gaussian.test,formula = ~G|A))
})

#delete file temporary file
if (file.exists("model.bug")) {
  file.remove("model.bug")}
