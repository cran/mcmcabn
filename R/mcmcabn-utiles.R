############################################################################### mcmcabn-utiles.R --- Author : Gilles Kratzer Document created : 13/02/2018 Last modification :

##-------------------------------------------------------------------------
## Internal function to test input of mcmcabn()
##-------------------------------------------------------------------------

.tests.mcmcabn <- function(score.cache, data.dists, max.parents, mcmc.scheme, seed, verbose, start.dag, prior.dag, prior.lambda,
    prob.rev, prob.mbr, prior.choice) {
    # start tests
    if (is.null(score.cache))
        stop("A cache of score should be provided. You can produce it using the R package abn.")

    if (max(rowSums(score.cache$node.defn)) > (max.parents+1))
        stop("Check max.parents. It should be the same as the one used in abn::buildscorecache() R function")

    if (length(mcmc.scheme) != 3)
        stop("An MCMC scheme have to be provided. It should be such that c(returned,thinned,burned) made of non negative integers.")

    if (!is.numeric(mcmc.scheme[1]) | !is.numeric(mcmc.scheme[2]) | !is.numeric(mcmc.scheme[3]) | mcmc.scheme[1] < 0 |
        mcmc.scheme[2] < 0 | mcmc.scheme[3] < 0) {
        stop("An MCMC scheme have to be provided. It should be such that c(returned,thinned,burned) made of non negative integers.")
    }

    if (max.parents < 1 || max.parents > length(data.dists))
        stop("max.parents makes no sense.")

    if (is.numeric(prob.rev) && prob.rev > 1 && prob.rev < 0)
        stop("prob.rev should be a probability.")

    if (is.numeric(prob.mbr) && prob.mbr > 1 && prob.mbr < 0)
        stop("prob.mbr should be a probability.")

    if (is.matrix(start.dag) && (dim(start.dag)[1] != length(data.dists)))
        stop("start.dag should be a squared matrix with dimension equal to the number of variables.")

    if (is.matrix(start.dag) && is.null(topoSortMAT((start.dag)))) {
        stop("start.dag should be a named DAG.")
    }

    if (is.matrix(start.dag) && max(rowSums(start.dag))>max.parents) {
      stop("start.dag should not have more parent than max.parent argument.")
    }

    if (length(data.dists) != max(score.cache$children))
        stop("data.dists should be a named list of all variables used to build the cache of precomputed score.")

    if (length(data.dists) == 0 || is.null(names(data.dists))) {
        stop("data.dists should be a named list of all variables used to build the cache of precomputed score.")
    }

    if (prior.choice %!in% 1:3) {
        stop("prior.choice should be either 1,2 or 3.")
    }

    if (is.matrix(prior.dag) && dim(prior.dag)[1] != length(data.dists))
        stop("prior.dag should be a squared matrix with dimension equal to the number of variables.")
}

##-------------------------------------------------------------------------
## Internal function that call multiple times strsplit() and remove space
##-------------------------------------------------------------------------

strsplits <- function(x, splits, ...) {
    for (split in splits) {
        x <- unlist(strsplit(x, split, ...))
    }
    x <- gsub(" ", "", x, fixed = TRUE)  #remove space
    return(x[!x == ""])  # Remove empty values
}

##-------------------------------------------------------------------------
## Internal function that produce a square matrix length(name) with {0,1} depending on f. f have to start with ~
## terms are entries of name terms are separated by + term1 | term2 indicates col(term1) row(term2) puts a 1 term1 |
## term2:term3: ... : is used as a sep . = all terms in name
##-------------------------------------------------------------------------

formula.mcmcabn <- function(f, name) {

    name_orignial <- name

    f <- as.character(f)

    ## tests for consistence ----------------------------------------------------------------------

    ## transformation name + or | or : or . or name to name_name
    if (sum((c("+", "|", ":", ".") %in% unlist(strsplit(name, split = c(""))))) != 0) {
        for (i in 1:length(name)) {
            if (sum(unlist(strsplit(name[i], split = c(""))) %in% c("+")) != 0) {
                f[[2]] <- gsub(name[i], gsub("+", "_", name[i], fixed = TRUE), f[[2]], fixed = TRUE)
                name[i] <- gsub("+", "_", name[i], fixed = TRUE)
            }
            if (sum(unlist(strsplit(name[i], split = c(""))) %in% c("|")) != 0) {
                f[[2]] <- gsub(name[i], gsub("|", "_", name[i], fixed = TRUE), f[[2]], fixed = TRUE)
                name[i] <- gsub("|", "_", name[i], fixed = TRUE)
            }
            if (sum(unlist(strsplit(name[i], split = c(""))) %in% c(":")) != 0) {
                f[[2]] <- gsub(name[i], gsub(":", "_", name[i], fixed = TRUE), f[[2]], fixed = TRUE)
                name[i] <- gsub(":", "_", name[i], fixed = TRUE)
            }
            if (sum(unlist(strsplit(name[i], split = c(""))) %in% c(".")) != 0) {
                f[[2]] <- gsub(name[i], gsub(".", "_", name[i], fixed = TRUE), f[[2]], fixed = TRUE)
                name[i] <- gsub(".", "_", name[i], fixed = TRUE)
            }
        }
    }

    ## collapse name
    name.c <- paste(name, collapse = ":")
    ## Split by terms
    f.p <- strsplit(x = f[[2]], split = "+", fixed = TRUE)

    ## nothing more than name variable in the dag formula
    tmp.test <- strsplits(x = f[[2]], splits = c("+", "|", ":", "."), fixed = TRUE)
    if (sum(!(tmp.test %in% name)) != 0) {
        stop("Formula contains some variables not used in the mcmcabn search")
    }
    ## End of tests for consistence ----------------------------------------------------------------

    ## creat the void matrix
    out <- matrix(data = 0, nrow = length(name), ncol = length(name))

    ## delete all spaces
    f.p <- gsub(" ", "", f.p[[1]], fixed = TRUE)

    ## replace '.' by all names
    f.p.completed <- gsub(".", name.c, f.p, fixed = TRUE)

    ## atomization of left term


    ## contruction of the output matrix
    for (i in 1:length(f.p)) {
        tmp <- f.p.completed[i]

        ## forget unique terms -> test for |
        if (grepl("|", tmp, fixed = TRUE)) {

            ## split wrt |
            tmp.p <- strsplit(x = tmp, split = "|", fixed = TRUE)

            ## test for multiple terms and contruction of the list first term
            if (grepl(":", tmp.p[[1]][1])) {
                tmp.p.p.1 <- strsplit(x = tmp.p[[1]][1], split = ":", fixed = TRUE)
            }
            if (!grepl(":", tmp.p[[1]][1])) {
                tmp.p.p.1 <- tmp.p[[1]][1]
            }

            ## test for multiple terms and contruction of the list second term
            if (grepl(":", tmp.p[[1]][2])) {
                tmp.p.p.2 <- strsplit(x = tmp.p[[1]][2], split = ":", fixed = TRUE)
            }
            if (!grepl(":", tmp.p[[1]][2])) {
                tmp.p.p.2 <- tmp.p[[1]][2]
            }

            ## loop over the
            for (j in 1:length(tmp.p.p.1[[1]])) {
                for (k in 1:length(tmp.p.p.2[[1]])) {
                  ## update of matrix
                  out[grep(tmp.p.p.1[[1]][j], name), grep(tmp.p.p.2[[1]][k], name)] <- 1

                }
            }
        }

    }

    ## avoid auto dependance
    diag(out) <- 0

    ## only 0 and 1
    out[out > 1] <- 1

    ## naming
    colnames(out) <- name_orignial
    rownames(out) <- name_orignial
    ## output
    return(out)
}


##-------------------------------------------------------------------------
## Ancestor function
##-------------------------------------------------------------------------

##-------------------------------------------------------------------------
## descendent function
##-------------------------------------------------------------------------

descendents <- function(nodes, dag) {
    if (!identical(topoSortMAT(amat = dag, index = FALSE), character(0))) {
        if (!is.integer0(which(x = dag[, nodes] == 1, arr.ind = TRUE)) && sum(diag(dag)) == 0) {
            if (length(nodes) == 1) {
                tmp <- unname(which(x = dag[, nodes] == 1, arr.ind = TRUE))
            } else {
                tmp <- t(unname(which(x = dag[, nodes] == 1, arr.ind = TRUE)))[1, ]
            }
            return(unique(unname(c(tmp, descendents(nodes = (tmp), dag = dag)))))
        } else {
            return(NULL)
        }
    }
    return("Not a DAG")
}

##-------------------------------------------------------------------------
## range function
##-------------------------------------------------------------------------

range01 <- function(x) {
    (x - min(x))/(max(x) - min(x))
}

##-------------------------------------------------------------------------
## not in function
##-------------------------------------------------------------------------

"%!in%" <- function(x, y) !(x %in% y)

##-------------------------------------------------------------------------
## is.integer0
##-------------------------------------------------------------------------

is.integer0 <- function(x) {
    is.integer(x) && length(x) == 0L
}

##-------------------------------------------------------------------------
## scoring DAGs
##-------------------------------------------------------------------------

score.dag <- function(dag,bsc.score,sc){
  n.var <- dim(dag)[1]
  score <- 0
  for (a in 1:n.var) {
    sc.tmp <- sc[bsc.score$children == a, ]
    score <- sum(min(sc.tmp[(apply(sc.tmp, 1, function(x) identical(unname(x[1:n.var]), unname(dag[a,])))), n.var + 1]), score)
  }
  return(score)
}
