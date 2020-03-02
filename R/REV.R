REV <- function(n.var, dag.tmp, retain, ban, max.parents, sc, score.cache, score, verbose, heating) {

    rejection <- 1
    A <- 0

    # stor number of edges

    n.edges <- sum(dag.tmp)

    # randomly select one

    if (sum(dag.tmp*(1-retain)*(1-t(ban))) != 0){

            #i->j
            selected.edge <- which(x = ((dag.tmp*(1-retain)*(1-t(ban))) == 1), arr.ind = TRUE)[sample(x = 1:sum(dag.tmp*(1-retain)*(1-t(ban))), size = 1), ,drop=FALSE]

            # store current parent set (j)
            parent.set <- dag.tmp[selected.edge[1], ,drop=FALSE]

            # remove parent set (i and j)
            dag.M.dot <- dag.tmp
            dag.M.dot[selected.edge[1],] <- 0
            dag.M.dot[selected.edge[2],] <- 0

            # store descendents j row i col

            descendents.M.dot.i <- descendents(nodes = unname(selected.edge[2]), dag = dag.M.dot)
            descendents.M.dot.j <- descendents(nodes = unname(selected.edge[1]), dag = dag.M.dot)

            # mark all parent sets of i score that do not include j

            sc.tmp <- sc[score.cache$children == selected.edge[2], ,drop=FALSE]
            sc.tmp <- sc.tmp[sc.tmp[, selected.edge[1]] == 1, ,drop=FALSE]


                #if (!is.null(descendents.M.dot.i)) {

                     sc.tmp <- sc.tmp[rowSums(sc.tmp[, descendents.M.dot.i,drop=FALSE] == 0) == length(descendents.M.dot.i), ,drop=FALSE]
                #}

                  new.parent.i <- sc.tmp[sample(x = 1:nrow(sc.tmp), size = 1, prob = range01(sc.tmp[, ncol(sc.tmp)] -
                    sum(sc.tmp[, ncol(sc.tmp)]))), ,drop=FALSE]
                  new.parent.i <- new.parent.i[-length(new.parent.i)]
                  # store partition function
                  #z.star.x.i.M.dot <- (logSumExp(sc.tmp[, ncol(sc.tmp)]))
                  #z.star.x.i.M.dot <- (sum(sc.tmp[, ncol(sc.tmp)]))

                # M direct sum
                dag.M.cross <- dag.M.dot
                dag.M.cross[selected.edge[2], ] <- new.parent.i

                # descendant of i
                descendents.M.cross.i <- descendents(nodes = unname(selected.edge[2]), dag = dag.M.cross)
                descendents.M.cross.j <- descendents(nodes = unname(selected.edge[1]), dag = dag.M.cross)
                #descendents.M.j <- descendents(nodes = unname(selected.edge[1]), dag = dag.tmp)

                # score node j not descendant of M
                sc.tmp <- sc[score.cache$children == selected.edge[1], ,drop=FALSE]

                # remove descendents of i
                #if (!is.null(descendents.M.cross.j)) {
                    #sc.tmp <- sc.tmp[rowSums(sc.tmp[, descendents.M.j,drop=FALSE] == 0) == length(descendents.M.j),, drop=FALSE]
                    sc.tmp <- sc.tmp[rowSums(sc.tmp[, descendents.M.cross.j,drop=FALSE] == 0) == length(descendents.M.cross.j),, drop=FALSE]
                 #}

                    new.parent.j <- sc.tmp[sample(x = 1:nrow(sc.tmp), size = 1, prob = range01(sc.tmp[, ncol(sc.tmp)] -
                      sum(sc.tmp[, ncol(sc.tmp)]))), ,drop=FALSE]
                    new.parent.j <- new.parent.j[-length(new.parent.j)]

                    # store partition function
                    #z.x.j.M.cross <- (logSumExp(sc.tmp[, ncol(sc.tmp)]))
                    #z.x.j.M.cross <- (sum(sc.tmp[, ncol(sc.tmp)]))

                  # M tilde

                  dag.M.tilde <- dag.M.cross
                  dag.M.tilde[selected.edge[1],] <- new.parent.j

                  n.edges.tilde <- sum(dag.M.tilde)

                  ############################## computing acceptance probability ##############################################

                  if(FALSE){
                  # score node j that do not include i
                  sc.tmp <- sc[score.cache$children == selected.edge[1], ,drop=FALSE]
                  sc.tmp <- sc.tmp[sc.tmp[, selected.edge[2]] == 1, ,drop=FALSE]
                  sc.tmp <- sc.tmp[rowSums(sc.tmp[, descendents.M.dot.j,drop=FALSE] == 0) == length(descendents.M.dot.j), ,drop=FALSE]

                  #z.star.x.j.M.dot <- (logSumExp(sc.tmp[, ncol(sc.tmp)]))

                  dag.M.tilde.cross <- dag.M.dot
                  dag.M.tilde.cross[selected.edge[1],] <- parent.set

                  # descendant
                  descendents.M.tilde.cross.i <- descendents(nodes = unname(selected.edge[2]), dag = dag.M.tilde.cross)

                  sc.tmp <- sc[score.cache$children == selected.edge[2], ,drop=FALSE]
                  if (!is.null(descendents.M.tilde.cross.i)) {
                       sc.tmp <- sc.tmp[rowSums(sc.tmp[, descendents.M.tilde.cross.i,drop=FALSE] == 0) == length(descendents.M.tilde.cross.i),
                        ,drop=FALSE]
                    }
                    z.x.i.M.tilde.cross <- (logSumExp(sc.tmp[, ncol(sc.tmp)]))
                  }

                  ############################## Acceptance probability

                  #score.A <- min(exp((z.star.x.i.M.dot / (z.star.x.j.M.dot) * (z.x.j.M.cross) / (z.x.i.M.tilde.cross))*n.edges/n.edges.tilde),1)
                  s.proposed <- score.dag(dag.M.tilde,score.cache,sc)
                  s.current <- score.dag(dag.tmp,score.cache,sc)

                  A <- min(exp(( s.proposed - s.current) * (n.edges/n.edges.tilde) ), 1)

                  #if(is.nan(score.A)){score.A <- 0}

                  #if((score.A)<0){score.A <- 0}

                  #A <- min(1, score.A)

                  #if (rbinom(n = 1, size = 1, prob = A) == 1) {
                  if (runif(1)<(exponent(A,heating))) {
                    rejection <- 0
                    dag.tmp <- dag.M.tilde
                    score <- score.dag(dag.M.tilde,score.cache,sc)
                  }
    }#eoif


    ############################## Return

    return(list(dag.tmp = dag.tmp, score = score, alpha = A, rejection = rejection))
}  #EOF
