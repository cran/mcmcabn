##############################################################################
## mcmcabn.R --- Author : Gilles Kratzer Document created: 28/11/2019 -

##-------------------------------------------------------------------------
## Coupled heated MCMC procedure
##-------------------------------------------------------------------------

CoupledHeatedmcmcabn <- function(score.cache = NULL, score = "mlik", data.dists = NULL, max.parents = 1, mcmc.scheme = c(100, 1000,
    1000), seed = 42, verbose = FALSE, start.dag = NULL, prior.dag = NULL, prior.lambda = NULL, prob.rev = 0.05, prob.mbr = 0.05, heating = 1, n.chains = 4,
    prior.choice = 2) {

    #################################################### Tests

    .tests.mcmcabn(score.cache, data.dists, max.parents, mcmc.scheme, seed, verbose, start.dag, prior.dag, prior.lambda,
        prob.rev, prob.mbr, prior.choice,heating)
  if(n.chains<2){stop("Coupled heated MCMC runs should have two chains at minimum.")}

    ## end of tests

    ## format
    if (is.character(score))
        score <- tolower(score)
    score <- c("bic", "aic", "mdl", "mlik")[pmatch(score, c("bic", "aic", "mdl", "mlik"))]
    if (score %!in% c("bic", "aic", "mdl", "mlik"))
        stop("A method should be provided and be one of: bic, aic, mdl or mlik.")
    if (is.matrix(prior.dag)) {
        prior.choice <- 3
    }
    if (is.matrix(prior.dag) && is.null(prior.lambda)) {
        prior.lambda <- 1
    }

    n.var <- length(data.dists)

    prob.mc3 <- 1 - (prob.rev + prob.mbr)

    ## output
    out.dags <- array(data = NA, dim = c(n.var, n.var, mcmc.scheme[1] + 1))
    out.scores <- array(data = NA, dim = c(mcmc.scheme[1] + 1))

    out.scores.coupled <- array(data = NA, dim = c(mcmc.scheme[1] + 1, n.chains))

    out.alpha <- array(data = NA, dim = c(mcmc.scheme[1] + 1))
    out.method <- array(data = NA, dim = c(mcmc.scheme[1] + 1))
    out.rejection <- array(data = NA, dim = c(mcmc.scheme[1] + 1))
    heating.para <- array(data = NA, dim = c(n.chains))

    ## seeding
    set.seed(seed = seed)

    ## structural restriction

    retain <- score.cache$dag.retained
    ban <- score.cache$dag.banned

    if(is.null(retain)){retain <- matrix(data = 0, nrow = n.var, ncol = n.var)}

    if(is.null(ban)){ban <- matrix(data = 0, nrow = n.var, ncol = n.var)}


    ## Initializing matrix
    dag.tmp <- matrix(data = 0, nrow = n.var, ncol = n.var)

    ## start zero matrix if(!is.null(start.dag)){ colnames(dag.tmp) <- rownames(dag.tmp)<-name.dag<-sample(1:n.var) }

    ## start random matrix
    if (is.null(start.dag)) {
        start.dag <- "random"
    }
    if (is.character(start.dag) && start.dag == "random") {
        vec.tmp <- c(rep(1, max.parents), rep(0, 2 * n.var - max.parents))
        for (lines in 1:(n.var - 1)) {
            dag.tmp[1 + lines, 1:lines] <- sample(vec.tmp)[1:lines]
        }
        colnames(dag.tmp) <- rownames(dag.tmp) <- name.dag <- sample(1:n.var)
        dag.tmp <- dag.tmp[order(name.dag), order(name.dag)]
        #check constraints
        dag.tmp <- dag.tmp + retain
        dag.tmp <- dag.tmp*(1-ban)
        ## run a series of checks on the DAG passed
        dag.tmp <- abs(dag.tmp)
        ## consistency checks
        diag(dag.tmp) <- 0
        dag.tmp[dag.tmp > 0] <- 1
    }

    if (is.character(start.dag) && start.dag == "hc"){

      start.dag <- searchHeuristic(score.cache = score.cache,
                       score =  score,
                       num.searches = 100,
                       max.steps = 500,
                       seed = seed,
                       verbose = verbose,
                       start.dag = NULL,
                       dag.retained = NULL,
                       dag.banned = NULL,
                       algo = "hc",
                       tabu.memory = 10,
                       temperature = 0.9)

      start.dag <- start.dag[["dags"]][[which.max(x = unlist(start.dag$scores))]]

    }

    if (is.matrix(start.dag)) {
        start.dag[start.dag != 0] <- 1
        diag(start.dag) <- 0
        dag.tmp <- start.dag
        colnames(dag.tmp) <- rownames(dag.tmp) <- 1:n.var
    }


    ## scores
    if (score %in% c("bic", "aic", "mdl")) {
        tmp <- score.cache$node.defn[, as.numeric(colnames(dag.tmp))]
        colnames(tmp) <- as.numeric(colnames(dag.tmp))

        sc <- cbind(tmp, -score.cache[[score]])

    }
    if (score == "mlik") {

        tmp <- score.cache$node.defn[, as.numeric(colnames(dag.tmp))]
        colnames(tmp) <- as.numeric(colnames(dag.tmp))

        sc <- cbind(tmp, score.cache[[score]])
    }



    ## scoring init
    score.init <- score.dag(dag.tmp,score.cache,sc)

    ## out
    out.dags[, , 1] <- dag.tmp

    out.scores[1] <- score <- score.init

    out.scores.coupled[1,] <- rep(score.init, n.chains)

    out.alpha[1] <- alpha <- 0

    out.method[1] <- "MC3"

    out.rejection[1] <- 0

    for (j in 1:n.chains) {
      heating.para[j] <- 1/(1+heating*(j-1))
    }


    ## start mcmc search
    if (verbose)
        cat("Start MCMC burn in \n")

    j <- 1

    tmp.dag <- array(data = NA, dim = c(n.var, n.var, n.chains))
    tmp.scores <- array(data = NA, dim = c(n.chains))
    tmp.alpha <- array(data = NA, dim = c(n.chains))
    tmp.rejection <- array(data = NA, dim = c(n.chains))
    tmp.method <- array(data = NA, dim = c(n.chains))

#init
    tmp.dag[,,] <- dag.tmp
    tmp.scores[] <- score.init

    while (j <= mcmc.scheme[3]) {

        ## method choice:

for (c in 1:n.chains) {

  heating <- heating.para[c]
  dag.tmp <- tmp.dag[,,c]
  score <- tmp.scores[c]


  #print(score)
        method.choice <- sample(x = c("MC3", "REV", "MBR"), size = 1, prob = c(prob.mc3, prob.rev, prob.mbr))

        switch(method.choice, MC3 = {
            out <- MC3(n.var, (dag.tmp), retain, ban, max.parents, sc, score.cache, score, prior.choice, prior.lambda, prior.dag,
                verbose, heating)
            tmp.dag[,,c] <- out$dag.tmp
            tmp.scores[c] <- out$score
        }, REV = {
            if (verbose) {
                print("REV move")
            }
            out <- REV(n.var, (dag.tmp),retain, ban, max.parents, sc, score.cache, score, verbose, heating)
            tmp.dag[,,c] <- out$dag.tmp
            tmp.scores[c] <- out$score
        }, MBR = {
            if (verbose) {
                print("MBR move")
            }
            out <- MBR(n.var, dag.tmp,retain, ban, max.parents, sc, score.cache, score, verbose, heating)
            tmp.dag[,,c] <- out$dag.tmp
            tmp.scores[c] <- out$score
        })
}#eof: multiple heated chains

        ## swap between two randomly chosen chains
        rand.chains <- sample(x = 1:n.chains,2)

        alpha <- min(exponent(exp(tmp.scores[rand.chains[1]] - tmp.scores[rand.chains[2]]),heating.para[rand.chains[1]]) * exponent(exp(tmp.scores[rand.chains[2]] - tmp.scores[rand.chains[1]]),heating.para[rand.chains[2]]),1)
        if(is.nan(alpha)){alpha <- 0}

if (runif(1)<(alpha)) {
  tmp.dag[,,rand.chains[1]] <- tmp.dag[,,rand.chains[2]]
  tmp.scores[rand.chains[1]] <- tmp.scores[rand.chains[2]]
}
        j <- j + 1
    }  #EOF Burn in

    out.scores[1] <- tmp.scores[1]

    out.scores.coupled[1,] <- tmp.scores

    out.dags[,,1] <- tmp.dag[,,1]

    ### EOF: Burn-in phase!

    if (verbose)
        cat("Start MCMC search \n")
    for (mcmc.search in 1:mcmc.scheme[1]) {
        if (verbose)
            cat("processing search ...", mcmc.search, "\n")

        ## start blind search
        j <- 0


        while (j <= mcmc.scheme[2]) {

          method.choice <- sample(x = c("MC3", "REV", "MBR"), size = 1, prob = c(prob.mc3, prob.rev, prob.mbr))
          tmp.method[] <- method.choice

          for (c in 1:n.chains) {

            heating <- heating.para[c]
            dag.tmp <- tmp.dag[,,c]
            score <- tmp.scores[c]



            switch(method.choice, MC3 = {
                out <- MC3(n.var, (dag.tmp), retain, ban, max.parents, sc, score.cache, score, prior.choice, prior.lambda, prior.dag,
                  verbose, heating)
                tmp.dag[,,c] <- out$dag.tmp
                tmp.scores[c] <- out$score
                #tmp.alpha[c] <- out$alpha
                tmp.rejection[c] <- out$rejection
            }, REV = {
                if (verbose) {
                  print("REV move")
                }
                out <- REV(n.var, (dag.tmp),retain, ban, max.parents, sc, score.cache, score, verbose, heating)
                tmp.dag[,,c] <- out$dag.tmp
                tmp.scores[c] <- out$score
                #tmp.alpha[c] <- out$alpha
                tmp.rejection[c] <- out$rejection
            }, MBR = {
                if (verbose) {
                  print("MBR move")
                }
                out <- MBR(n.var, (dag.tmp), retain, ban, max.parents, sc, score.cache, score, verbose, heating)
                tmp.dag[,,c] <- out$dag.tmp
                tmp.scores[c] <- out$score
                #tmp.alpha[c] <- out$alpha
                tmp.rejection[c] <- out$rejection
            })
          }#eof: multiple heated chains

            ## swap between two randomly chosen chains
            rand.chains <- sample(x = 1:n.chains,2)

            alpha <- min(exponent(exp(tmp.scores[rand.chains[1]] - tmp.scores[rand.chains[2]]),heating.para[rand.chains[1]]) * exponent(exp(tmp.scores[rand.chains[2]] - tmp.scores[rand.chains[1]]),heating.para[rand.chains[2]]),1)
            if(is.nan(alpha)){alpha <- 0}


                        if (runif(1)<(alpha)) {
              storage.tmp1 <- tmp.dag[,,rand.chains[1]]
              storage.tmp2 <-tmp.scores[rand.chains[1]]
              storage.tmp3 <- tmp.method[rand.chains[1]]

              tmp.dag[,,rand.chains[1]] <- tmp.dag[,,rand.chains[2]]
              tmp.scores[rand.chains[1]] <- tmp.scores[rand.chains[2]]
              tmp.method[rand.chains[1]] <- tmp.method[rand.chains[2]]

               tmp.dag[,,rand.chains[2]] <- storage.tmp1
               tmp.scores[rand.chains[2]] <- storage.tmp2
               tmp.method[rand.chains[2]] <- storage.tmp3

              #rejection <- 0
            }

          j <- j + 1
        }  #EOF mcmc search


        out.dags[, , mcmc.search + 1] <- tmp.dag[,,1]

        out.scores[mcmc.search + 1] <- tmp.scores[1]

        out.alpha[mcmc.search + 1] <- alpha

        out.method[mcmc.search + 1] <- tmp.method[1]

        out.rejection[mcmc.search + 1] <- tmp.rejection[1]

        out.scores.coupled[mcmc.search + 1,] <- tmp.scores

        out.scores.coupled <- data.frame(out.scores.coupled)

        names(out.scores.coupled) <- 1:n.chains

    }  #EOF mcmc search
    out <- list(dags = out.dags, scores = out.scores, alpha = out.alpha, method = out.method, rejection = out.rejection,
        iterations = mcmc.scheme[1] * (mcmc.scheme[2] + 1), thinning = mcmc.scheme[2], burnin = mcmc.scheme[3], data.dist = data.dists, heating = heating.para, score.coupled = out.scores.coupled)

    class(out) <- "mcmcabn"

    return(out)

}  #eof
