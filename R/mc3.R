mc3 <- function(n.var, dag.tmp, retain, ban, max.parents, sc, score.cache, score, prior.choice, prior.lambda, prior.dag, verbose, heating) {

    ## construction of neighbours list
    neighbours.list <- NULL
    alpha <- 0
    rejection <- 1

    ## choose mc move
    if (sum(dag.tmp) > 0) {
        if (sum(dag.tmp) < n.var * max.parents) {
            mc.move <- sample(x = c("removal", "addition", "reversal"), size = 1)
        } else {
            mc.move <- "removal"  #rep
        }
    } else {
        mc.move <- "addition"  #rep
    }

    switch(mc.move, removal = {
        ## Removal of an arc
        for (a in 1:n.var) {
            for (b in 1:n.var) {
                if (dag.tmp[a, b] == 1 && retain[a, b]!=1) {
                  tmp <- dag.tmp
                  tmp[a, b] <- 0
                  neighbours.list[[length(neighbours.list) + 1]] <- tmp
                }
            }
        }
    }, addition = {

        ## Addition of an arc
        for (a in 1:n.var) {
            for (b in 1:n.var) {
                if(((dag.tmp[a, b] == 0 & dag.tmp[b, a] == 0) & sum(dag.tmp[a, ]) <= (max.parents - 1)) && ban[a,b]!=1) {
                  tmp <- dag.tmp
                  # print(tmp)
                  tmp[a, b] <- 1


                  if (!identical(topoSortMAT(amat = tmp, index = FALSE), character(0))) {
                    neighbours.list[[length(neighbours.list) + 1]] <- tmp
                  }
                }
            }
        }
    }, reversal = {

        ## Reversal of an arc
        for (a in 1:n.var) {
            for (b in 1:n.var) {
                if ((dag.tmp[a, b] == 1 & sum(dag.tmp[b, ]) <= (max.parents - 1)) && retain[a,b] && ban[b,a]!=1) {
                  tmp <- dag.tmp
                  tmp[a, b] <- 0
                  tmp[b, a] <- 1
                  if (!identical(topoSortMAT(amat = tmp, index = FALSE), character(0))) {
                    neighbours.list[[length(neighbours.list) + 1]] <- tmp
                  }
                }
            }
        }
    })


    ## choose of an mcmc move
    n.G <- length(neighbours.list)
    if (0 != n.G) {
        mcmc.move <- sample(x = 1:n.G, size = 1, replace = FALSE)
        dag.gprime <- neighbours.list[[mcmc.move]]

        ## construction of neighbours list of GPRIME
        neighbours.list.gprime <- NULL

        switch(mc.move, removal = {
            ## Removal of an arc
            for (a in 1:n.var) {
                for (b in 1:n.var) {
                  if (dag.gprime[a, b] == 1) {
                    tmp <- dag.gprime
                    tmp[a, b] <- 0
                    neighbours.list.gprime[[length(neighbours.list.gprime) + 1]] <- tmp
                  }
                }
            }
        }, addition = {

            ## Addition of an arc
            for (a in 1:n.var) {
                for (b in 1:n.var) {
                  if ((dag.tmp[a, b] == 0 & dag.tmp[b, a] == 0) & sum(dag.tmp[a, ]) <= (max.parents - 1)) {
                    tmp <- dag.gprime
                    # print(tmp)
                    tmp[a, b] <- 1
                    if (!identical(topoSortMAT(amat = tmp, index = FALSE), character(0))) {
                      neighbours.list.gprime[[length(neighbours.list.gprime) + 1]] <- tmp
                    }
                  }
                }
            }
        }, reversal = {

            ## Reveral of an arc
            for (a in 1:n.var) {
                for (b in 1:n.var) {
                  if (dag.tmp[a, b] == 1 & sum(dag.tmp[b, ]) <= (max.parents - 1)) {
                    tmp <- dag.gprime
                    tmp[a, b] <- 0
                    tmp[b, a] <- 1
                    if (!identical(topoSortMAT(amat = tmp, index = FALSE), character(0))) {
                      neighbours.list.gprime[[length(neighbours.list.gprime) + 1]] <- tmp
                    }
                  }
                }
            }
        })

        n.Gprime <- length(neighbours.list.gprime)

        ## prior computation
        prior.G <- prior.Gprime <- 1
        score.G <- score.Gprime <- 0

        # user prior
        if (prior.choice == 3) {
            prior.G <- exp(-prior.lambda * sum(abs(dag.tmp - prior.dag)))
            prior.Gprime <- exp(-prior.lambda * sum(abs(dag.gprime - prior.dag)))
        }

        for (a in 1:n.var) {
            # Koivisto priors
            if (prior.choice == 2) {
                prior.G <- prod(1/choose(n = (n.var - 1), k = sum(dag.tmp[a, ])), prior.G)
                prior.Gprime <- prod(1/choose(n = (n.var - 1), k = sum(dag.gprime[a, ])), prior.Gprime)
            }

            sc.tmp <- sc[score.cache$children == a,,drop = FALSE]
            score.G <- sum(min(sc.tmp[(apply(sc.tmp, 1, function(x) identical(unname(x[1:n.var]), unname(dag.tmp[a,
                ])))), n.var + 1]), score.G)
            score.Gprime <- sum(min(sc.tmp[(apply(sc.tmp, 1, function(x) identical(unname(x[1:n.var]), unname(dag.gprime[a,
                ])))), n.var + 1]), score.Gprime)
        }


        score.Gprime.scaled <- score.dag(dag.gprime,score.cache,sc)
        score.G.scaled <- score.dag(dag.tmp,score.cache,sc)

        alpha <- min(exp( (score.Gprime.scaled - score.G.scaled) * (n.G/n.Gprime) * (prior.Gprime/prior.G)), 1)

        if(!is.numeric(alpha) | is.nan(alpha)) alpha <- 0

        score <- score.G

        if (!is.null(dag.gprime) && runif(1)<exponent(alpha,heating)) {
            dag.tmp <- dag.gprime
            score <- score.Gprime
            rejection <- 0
        }
    }

    return(list(dag.tmp = dag.tmp, score = score, alpha = alpha, rejection = rejection))
}  #EOF
