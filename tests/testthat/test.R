###############################################################################
## test.R ---
## Author          : Gilles Kratzer
## Document created: 25/02/2019
###############################################################################

##Purpose: Test the mcmcabn software

Sys.setenv("R_TESTS" = "")
#library(testthat)
#library(abn)
data(asia, package='bnlearn')
data(marks, package='bnlearn')
data(gaussian.test, package='bnlearn')

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

  out.sim.0 <- simulateAbn(data.dists = dist,n.chains = 1,n.adapt = 1000,n.thin = 1,n.iter = 1000,data.param = 0.4*data.param, simulate = TRUE,seed = 132,verbose = FALSE)

  bsc.compute.0 <- buildScoreCache(data.df = out.sim.0, data.dists = dist,  max.parents = 3)

  mc3.out <- mcmcabn(score.cache = bsc.compute.0,
                  score = "mlik",
                  data.dists = dist,
                  max.parents = 3,
                  mcmc.scheme = c(50,0,0),
                  seed = 42,
                  verbose = FALSE,
                  start.dag = "random",
                  prob.rev = 0,
                  prob.mbr = 0,
                  prior.choice = 1,heating = 0.5)

  dag <- mostProbable(score.cache = bsc.compute.0,prior.choice = 1)
  expect_equal(max(mc3.out$scores),fitAbn(object = dag,data.df = out.sim.0,data.dists = dist)$mlik)

  mc3.out <- mcmcabn(score.cache = bsc.compute.0,
                     score = "mlik",
                     data.dists = dist,
                     max.parents = 3,
                     mcmc.scheme = c(50,0,0),
                     seed = 465,
                     verbose = FALSE,
                     start.dag = "random",
                     prob.rev = 0,
                     prob.mbr = 0,
                     prior.choice = 2)

  dag <- mostProbable(score.cache = bsc.compute.0,prior.choice = 2)
  expect_equal(max(mc3.out$scores),fitAbn(object = dag,data.df = out.sim.0,data.dists = dist)$mlik, tol=0.0001)

  #test influence of user define prior
  data.param.eq <- matrix(data = 0,nrow = 5,ncol = 5)
  mc3.out <- mcmcabn(score.cache = bsc.compute.0,
                     score = "mlik",
                     data.dists = dist,
                     max.parents = 3,
                     mcmc.scheme = c(100,0,0),
                     seed = 165654,
                     verbose = FALSE,
                     start.dag = "random",
                     prob.rev = 0,
                     prob.mbr = 0,
                     prior.dag = data.param.eq,
                     prior.lambda = 10000,
                     prior.choice = 3)

  expect_false(table(apply(mc3.out$dags,3, sum))[1]<50)

  data.param.eq <- matrix(data = 1,nrow = 5,ncol = 5)
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

  out.sim.0 <- invisible(simulateAbn(data.dists = dist,n.chains = 1,n.adapt = 1000,n.thin = 1,n.iter = 1000,data.param = 0.4*data.param, simulate = TRUE,seed = 132))

  bsc.compute.0 <- buildScoreCache(data.df = out.sim.0, data.dists = dist,  max.parents = 3)

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

  dag <- mostProbable(score.cache = bsc.compute.0,prior.choice = 1)
  expect_equal(max(mc3.out$scores),fitAbn(object = dag,data.df = out.sim.0,data.dists = dist)$mlik)

  expect_silent(a<-mcmcabn(score.cache = bsc.compute.0,
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

  ##test theoretical

  dist <- list(A="binomial",
               B="binomial",
               C="binomial")

  data.param <- matrix(data = c(0,1,0,0,0,0,0,0,0),nrow = 3L,ncol = 3L,byrow = TRUE)

  diag(data.param) <- 1
  colnames(data.param) <- rownames(data.param) <- names(dist)

  out.sim.0 <- invisible(simulateAbn(data.dists = dist,n.chains = 1,n.adapt = 1000,n.thin = 1,n.iter = 10,data.param = 0.4*data.param, simulate = TRUE,seed = 132))

  bsc.compute.0 <- buildScoreCache(data.df = out.sim.0, data.dists = dist,  max.parents = 3)

  #
  start.m <- matrix(data = c(0,0,0,1,0,0,0,0,0),nrow = 3L,ncol = 3L,byrow = TRUE)

  retain.m <- matrix(data = c(0,0,0,0,0,0,0,0,0),nrow = 3L,ncol = 3L,byrow = TRUE)

  ban.m <- matrix(data = c(0,0,0,0,0,0,0,0,0),nrow = 3L,ncol = 3L,byrow = TRUE)

  colnames(start.m) <- rownames(start.m) <- names(dist)

  tmp <- bsc.compute.0$node.defn[, 1:3]
  colnames(tmp) <- 1:3
  sc <- cbind(tmp, bsc.compute.0[["mlik"]])

  mcmcabn:::REV(n.var = 3L,dag.tmp = start.m,retain = retain.m,ban = ban.m,max.parents = 3,sc = sc,score.cache = bsc.compute.0,score = -10,verbose = TRUE,heating = 1)

  start.m <- matrix(data = c(0,1,0,1,0,1,0,0,0),nrow = 3L,ncol = 3L,byrow = TRUE)

  mcmcabn:::REV(n.var = 3L,dag.tmp = start.m,retain = retain.m,ban = ban.m,max.parents = 3,sc = sc,score.cache = bsc.compute.0,score = -10,verbose = TRUE,heating = 1)
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

  out.sim.0 <- invisible(simulateAbn(data.dists = dist,n.chains = 1,n.adapt = 1000,n.thin = 1,n.iter = 1000,data.param = 0.4*data.param, simulate = TRUE,seed = 132))

  bsc.compute.0 <- buildScoreCache(data.df = out.sim.0, data.dists = dist,  max.parents = 3)

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

  dag <- mostProbable(score.cache = bsc.compute.0,prior.choice = 1)
  expect_equal(max(mc3.out$scores),fitAbn(object = dag,data.df = out.sim.0,data.dists = dist)$mlik)

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

  out.sim.0 <- invisible(simulateAbn(data.dists = dist,n.chains = 1,n.adapt = 20,n.thin = 1,n.iter = 100,data.param = data.param.0, simulate = TRUE,seed = 132))

  bsc.compute.0 <- buildScoreCache(data.df = out.sim.0, data.dists = dist,  max.parents = 2)

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
  #data(asia)

  dist.asia <- list(Asia = "binomial",
                    Smoking = "binomial",
                    Tuberculosis = "binomial",
                    LungCancer = "binomial",
                    Bronchitis = "binomial",
                    Either = "binomial",
                    XRay = "binomial",
                    Dyspnea = "binomial")

  colnames(asia) <- c("Asia","Smoking", "Tuberculosis", "LungCancer", "Bronchitis", "Either", "XRay", "Dyspnea")


  bsc.compute.asia <- buildScoreCache(data.df = asia,
                                      data.dists = dist.asia,
                                      max.parents = 2)


  mcmc.out.asia <- mcmcabn(score.cache = bsc.compute.asia,
                           score = "mlik",
                           data.dists = dist.asia,
                           max.parents = 2,
                           mcmc.scheme = c(500,0,0),
                           seed = 45235,
                           verbose = FALSE,
                           start.dag = "random",
                           prob.rev = 0.2,
                           prob.mbr = 0.2,
                           prior.choice = 1,heating = 0.7)


  #maximum scoring network using exact search (not MCMC based)
  dag <- mostProbable(score.cache = bsc.compute.asia)
  expect_equal(max(mcmc.out.asia$scores),fitAbn(object = dag,data.df = asia,data.dists = dist.asia)$mlik)

  mcmc.out.asia <- mcmcabn(score.cache = bsc.compute.asia,
                           score = "mlik",
                           data.dists = dist.asia,
                           max.parents = 2,
                           mcmc.scheme = c(150,0,0),
                           seed = 124019,
                           verbose = FALSE,
                           start.dag = "hc",
                           prob.rev = 0.2,
                           prob.mbr = 0.2,
                           prior.choice = 1,heating = 0.7)

  expect_equal(max(mcmc.out.asia$scores),fitAbn(object = dag,data.df = asia,data.dists = dist.asia)$mlik)

  ## marks datasets

  #data(marks)

  dist.marks <- list(MECH = "gaussian",
                    VECT = "gaussian",
                    ALG = "gaussian",
                    ANL = "gaussian",
                    STAT = "gaussian")

  #colnames(asia) <- c("Asia","Smoking", "Tuberculosis", "LungCancer", "Bronchitis", "Either", "XRay", "Dyspnea")


  bsc.compute.marks <- buildScoreCache(data.df = marks,
                                      data.dists = dist.marks,
                                      max.parents = 2)


  mcmc.out.marks <- mcmcabn(score.cache = bsc.compute.marks,
                           score = "mlik",
                           data.dists = dist.marks,
                           max.parents = 2,
                           mcmc.scheme = c(150,0,0),
                           seed = 789,
                           verbose = FALSE,
                           start.dag = "random",
                           prob.rev = 0.1,
                           prob.mbr = 0.1,
                           prior.choice = 2)


  #maximum scoring network using exact search (not MCMC based)
  dag <- mostProbable(score.cache = bsc.compute.marks)
  expect_equal(max(mcmc.out.marks$scores),fitAbn(object = dag,data.df = marks,data.dists = dist.marks)$mlik)

  ##tests
  #data(gaussian.test)

  dist.gaussian.test <- list(A = "gaussian",
                     B = "gaussian",
                     C = "gaussian",
                     D = "gaussian",
                     E = "gaussian",
                     G = "gaussian",
                     H ="gaussian")

  colnames(gaussian.test) <- c("A","B","C","D","E","G","H")

  bsc.compute.gaussian.test <- buildScoreCache(data.df = gaussian.test,
                                       data.dists = dist.gaussian.test,
                                       max.parents = 2)


  mcmc.out.gaussian.test <- mcmcabn(score.cache = bsc.compute.gaussian.test,
                            score = "mlik",
                            data.dists = dist.gaussian.test,
                            max.parents = 2,
                            mcmc.scheme = c(1500,0,0),
                            seed = 148695,
                            verbose = FALSE,
                            start.dag = "random",
                            prob.rev = 0.2,
                            prob.mbr = 0.2,
                            prior.choice = 2)



  #maximum scoring network using exact search (not MCMC based)
  dag <- mostProbable(score.cache = bsc.compute.gaussian.test)
  expect_equal(max(mcmc.out.gaussian.test$scores),fitAbn(object = dag,data.df = gaussian.test,data.dists = dist.gaussian.test)$mlik, tol=10)

  })

test_that("query",{

  #data(gaussian.test)

  dist.gaussian.test <- list(A = "gaussian",
                             B = "gaussian",
                             C = "gaussian",
                             D = "gaussian",
                             E = "gaussian",
                             G = "gaussian",
                             H ="gaussian")

  colnames(gaussian.test) <- c("A","B","C","D","E","G","H")

  bsc.compute.gaussian.test <- buildScoreCache(data.df = gaussian.test,
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

test_that("mcmcabnHeated",{
#data(asia)

dist.asia <- list(Asia = "binomial",
                  Smoking = "binomial",
                  Tuberculosis = "binomial",
                  LungCancer = "binomial",
                  Bronchitis = "binomial",
                  Either = "binomial",
                  XRay = "binomial",
                  Dyspnea = "binomial")

colnames(asia) <- c("Asia","Smoking", "Tuberculosis", "LungCancer", "Bronchitis", "Either", "XRay", "Dyspnea")


bsc.compute.asia <- buildScoreCache(data.df = asia,
                                    data.dists = dist.asia,
                                    max.parents = 2)


mcmc.out.asia <- CoupledHeatedmcmcabn(score.cache = bsc.compute.asia,
                         score = "mlik",
                         data.dists = dist.asia,
                         max.parents = 2,
                         mcmc.scheme = c(2000,0,0),
                         seed = 9,
                         verbose = FALSE,
                         start.dag = "random",
                         prob.rev = 0.1,
                         prob.mbr = 0.1,
                         prior.choice = 1,heating = 5,n.chains = 4)

#maximum scoring network using exact search (not MCMC based)
dag <- mostProbable(score.cache = bsc.compute.asia)
expect_equal(max(mcmc.out.asia$score.coupled),fitAbn(object = dag,data.df = asia,data.dists = dist.asia)$mlik, tol=5)


dist.marks <- list(MECH = "gaussian",
                   VECT = "gaussian",
                   ALG = "gaussian",
                   ANL = "gaussian",
                   STAT = "gaussian")

#colnames(asia) <- c("Asia","Smoking", "Tuberculosis", "LungCancer", "Bronchitis", "Either", "XRay", "Dyspnea")


bsc.compute.marks <- buildScoreCache(data.df = marks,
                                     data.dists = dist.marks,
                                     max.parents = 2)


mcmc.out.marks <- CoupledHeatedmcmcabn(score.cache = bsc.compute.marks,
                          score = "mlik",
                          data.dists = dist.marks,
                          max.parents = 2,
                          mcmc.scheme = c(1000,0,0),
                          seed = 415654,
                          verbose = FALSE,
                          start.dag = "random",
                          prob.rev = 0.1,
                          prob.mbr = 0.1,n.chains = 4,heating = 5,
                          prior.choice = 1)

#maximum scoring network using exact search (not MCMC based)
dag <- mostProbable(score.cache = bsc.compute.marks)
expect_equal(max(mcmc.out.marks$scores),fitAbn(object = dag,data.df = marks,data.dists = dist.marks)$mlik, tol = 0.5)

})

#delete file temporary file
if (file.exists("model.bug")) {
  file.remove("model.bug")}
