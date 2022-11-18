################################################################################ mcmcabn.R --- Author : Gilles Kratzer Document created: 22/10/2019 -

##-------------------------------------------------------------------------
## MCMC procedure
##-------------------------------------------------------------------------

mcmcabn <- function(score.cache = NULL, score = "mlik", data.dists = NULL, max.parents = 1, mcmc.scheme = c(100, 1000,
    1000), seed = 42, verbose = FALSE, start.dag = NULL, prior.dag = NULL, prior.lambda = NULL, prob.rev = 0.05, prob.mbr = 0.05, heating = 1,
    prior.choice = 2) {

    #################################################### Tests

    .tests.mcmcabn(score.cache, data.dists, max.parents, mcmc.scheme, seed, verbose, start.dag, prior.dag, prior.lambda,
        prob.rev, prob.mbr, prior.choice,heating)

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

    heating.orig <- NULL

    ## output
    out.dags <- array(data = NA, dim = c(n.var, n.var, mcmc.scheme[1] + 1))
    out.scores <- rep(NA, length = mcmc.scheme[1] + 1)
    out.alpha <- rep(NA, length = mcmc.scheme[1] + 1)
    out.method <- rep(NA, length = mcmc.scheme[1] + 1)
    out.rejection <- rep(NA, length = mcmc.scheme[1] + 1)
    out.heating <- rep(NA, length = mcmc.scheme[1] + 1)

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

    ## heating scheme
    if(is.numeric(heating) & heating>1){
      heating.orig <- heating
      heating <- 1/heating.orig
    }


    ## out
    out.dags[, , 1] <- dag.tmp

    out.scores[1] <- score <- score.init

    out.alpha[1] <- alpha <- 0

    out.method[1] <- "MC3"

    out.rejection[1] <- 0

    out.heating[1] <- heating

    ## start mcmc search
    if (verbose)
        cat("Start MCMC burn in \n")

    j <- 1


    while (j <= mcmc.scheme[3]) {

        ## method choice:

        method.choice <- sample(x = c("MC3", "REV", "MBR"), size = 1, prob = c(prob.mc3, prob.rev, prob.mbr))

        if(heating<1 & !is.null(heating.orig)){heating <-1/(heating.orig-j+1)}


        switch(method.choice, MC3 = {
            out <- MC3(n.var, (dag.tmp), retain, ban, max.parents, sc, score.cache, score, prior.choice, prior.lambda, prior.dag,
                verbose, heating)
            dag.tmp <- out$dag.tmp
            score <- out$score
        }, REV = {
            if (verbose) {
                print("REV move")
            }
            out <- REV(n.var, (dag.tmp),retain, ban, max.parents, sc, score.cache, score, verbose, heating)
            dag.tmp <- out$dag.tmp
            score <- out$score
        }, MBR = {
            if (verbose) {
                print("MBR move")
            }
            out <- MBR(n.var, dag.tmp,retain, ban, max.parents, sc, score.cache, score, verbose, heating)
            dag.tmp <- out$dag.tmp
            score <- out$score
        })

        j <- j + 1
    }  #EOF Burn in

    out.scores[1] <- score
    out.dags[, , 1] <- dag.tmp

    ### EOF: Burn-in phase!

    if(!is.null(heating.orig) & heating<1){
      heating.orig <- round(1/heating)
    }

    if (verbose)
        cat("Start MCMC search \n")
    for (mcmc.search in 1:mcmc.scheme[1]) {
        if (verbose)
            cat("processing search ...", mcmc.search, "\n")

        ## start blind search
        j <- 0

        if(heating<1 & !is.null(heating.orig)){heating <- 1/(heating.orig-mcmc.search+1)}

        while (j <= mcmc.scheme[2]) {

            method.choice <- sample(x = c("MC3", "REV", "MBR"), size = 1, prob = c(prob.mc3, prob.rev, prob.mbr))

            switch(method.choice, MC3 = {
                out <- MC3(n.var, (dag.tmp), retain, ban, max.parents, sc, score.cache, score, prior.choice, prior.lambda, prior.dag,
                  verbose, heating)
                dag.tmp <- out$dag.tmp
                score <- out$score
                alpha <- out$alpha
                rejection <- out$rejection
            }, REV = {
                if (verbose) {
                  print("REV move")
                }
                out <- REV(n.var, (dag.tmp),retain, ban, max.parents, sc, score.cache, score, verbose, heating)
                dag.tmp <- out$dag.tmp
                score <- out$score
                alpha <- out$alpha
                rejection <- out$rejection
            }, MBR = {
                if (verbose) {
                  print("MBR move")
                }
                out <- MBR(n.var, (dag.tmp), retain, ban, max.parents, sc, score.cache, score, verbose, heating)
                dag.tmp <- out$dag.tmp
                score <- out$score
                alpha <- out$alpha
                rejection <- out$rejection
            })

            j <- j + 1
        }  #EOF mcmc search

        out.dags[, , mcmc.search + 1] <- dag.tmp

        out.scores[mcmc.search + 1] <- score

        out.alpha[mcmc.search + 1] <- alpha

        out.method[mcmc.search + 1] <- method.choice

        out.rejection[mcmc.search + 1] <- rejection

        out.heating[mcmc.search + 1] <- heating

    }  #EOF mcmc search
    out <- list(dags = out.dags, scores = out.scores, alpha = out.alpha, method = out.method, rejection = out.rejection,
        iterations = mcmc.scheme[1] * (mcmc.scheme[2] + 1), thinning = mcmc.scheme[2], burnin = mcmc.scheme[3], data.dist = data.dists, heating = out.heating)

    class(out) <- "mcmcabn"

    return(out)

}  #eof
