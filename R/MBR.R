MBR <- function(n.var, dag.tmp, retain, ban, max.parents, sc, score.cache, score, verbose, heating) {

    ## scores
    score.G <- 0
    score.G.prime <- 0
    rejection <- 1
    A <- NULL

    # store number of edges

    n.edges <- sum(dag.tmp)

    # Randomly select a node Xi in current graph G, withprobability N^−1

    selected.node <- sample(x = 1:n.var, size = 1)

    # store current parent set

    parent.set <- dag.tmp[selected.node, ,drop=FALSE]

    # remove parent set + co parent
    dag.G.0 <- dag.tmp

    # descendant of node i

    descendents.i <- descendents(nodes = selected.node, dag = dag.G.0)

    #delete parent of children of node i
    dag.G.0[descendents.i, -selected.node] <- 0

    dag.G.0[selected.node, ] <- 0

    # from score take out parent + descendant
    sc.tmp <- sc[score.cache$children == selected.node, ,drop=FALSE]

    sc.tmp <- sc.tmp[rowSums(sc.tmp[, c(as.logical(parent.set), FALSE),drop=FALSE]) == sum(parent.set), ,drop=FALSE]
    sc.tmp <- sc.tmp[rowSums(sc.tmp[, descendents.i, drop=FALSE] == 0) == length(descendents.i), ,drop=FALSE]

    # sample new parent set

    if(nrow(sc.tmp)==0){new.parent.i <- matrix(data = 0,nrow =  1,ncol = n.var+1)}
    if(nrow(sc.tmp)==1){new.parent.i <- sc.tmp
    }else{
        new.parent.i <- sc.tmp[sample(x = 1:nrow(sc.tmp), size = 1, prob = range01(sc.tmp[, ncol(sc.tmp)] - sum(sc.tmp[, ncol(sc.tmp)]))), ,drop=FALSE]
        }

    new.parent.i <- new.parent.i[,-length(new.parent.i)]


    #z.i.G.0 <- (logSumExp(sc.tmp[, ncol(sc.tmp)]))


            dag.G.1 <- dag.G.0
            dag.G.1[selected.node, ] <- new.parent.i

            ## new co parent set

            if (sum(dag.G.1[, selected.node]) != 0) {

                score.G <- numeric(sum(dag.G.1[, selected.node]))
                # randomization of children
                order.child <- sample(1:sum(dag.G.1[, selected.node]))
                for (j in order.child) {
                  # descendant of node j

                  descendents.j <- descendents(nodes = j, dag = dag.G.1)


                    # from score take out parent + descendant
                    sc.tmp <- sc[score.cache$children == j, , drop=FALSE]
                    sc.tmp <- sc.tmp[sc.tmp[, selected.node] == 1,,drop=FALSE]

                    sc.tmp <- sc.tmp[rowSums(sc.tmp[,descendents.j,drop=FALSE] == 0) == length(descendents.j),, drop=FALSE]

                    # sample a new parent set

                    if(nrow(sc.tmp)==0){new.parent.j <- matrix(data = 0,nrow =  1,ncol = n.var+1)}
                    if(nrow(sc.tmp)==1){new.parent.i <- sc.tmp}
                    if(nrow(sc.tmp)>1){
                        new.parent.j <- sc.tmp[sample(x = 1:nrow(sc.tmp), size = 1, prob = range01(sc.tmp[, ncol(sc.tmp)] - sum(sc.tmp[, ncol(sc.tmp)]))), ,drop=FALSE]
                        }

                    new.parent.j <- new.parent.j[-length(new.parent.j)]
                    score.G[j] <- sum(sc.tmp[, ncol(sc.tmp)])
                    #score.G[j] <- logSumExp(sc.tmp[, ncol(sc.tmp)])

                }
            }

            dag.MBR <- dag.G.1

            ############################ complementary inverse move
if(FALSE){

            # parent set
            parent.set.G.1 <- dag.G.1[selected.node, ]

            # delete co-parent

            tmp <- which(x = dag.G.1[, selected.node] == 1, arr.ind = TRUE)
            if (!is.integer0(tmp)) {
                dag.G.1[tmp[1], -selected.node] <- 0
            }

            # delete parents
            dag.G.1[selected.node, ] <- 0


            # from score take out parent + descendant
            sc.tmp <- sc[score.cache$children == selected.node, ,drop=FALSE]
            if (sum(parent.set.G.1) != 0) {

                  sc.tmp <- sc.tmp[rowSums(sc.tmp[, c(as.logical(parent.set.G.1), FALSE),drop=FALSE]) == sum(parent.set.G.1), ,drop=FALSE]

                # remove descendents of i
                if (!is.null(descendents.i)) {
                  if (length(descendents.i) == 1) {
                    if (is.null(dim(sc.tmp))) {
                      sc.tmp <- sc.tmp[sc.tmp[descendents.i] == 0]
                    } else {
                      sc.tmp <- sc.tmp[sc.tmp[, descendents.i] == 0, ,drop=FALSE]
                    }

                  } else {

                    if (is.null(dim(sc.tmp))) {
                      sc.tmp <- sc.tmp[sum(sc.tmp[descendents.i] == 0) == length(descendents.i)]
                    } else {
                      sc.tmp <- sc.tmp[rowSums(sc.tmp[, descendents.i,drop=FALSE] == 0) == length(descendents.i), ,drop=FALSE]
                    }
                  }
                }

                # store score
                if (is.null(ncol(sc.tmp))) {
                  z.star.i.G.prime.0 <- sc.tmp[length(sc.tmp)]

                } else {
                  z.star.i.G.prime.0 <- logSumExp(sc.tmp[, ncol(sc.tmp)])
                }

                # equivalent to G1
                dag.G.1[selected.node, ] <- parent.set

                if (sum(dag.G.1[, selected.node]) != 0) {

                  score.G.prime <- numeric(sum(dag.G.1[, selected.node]))

                  for (j in 1:sum(dag.G.1[, selected.node])) {
                    # descendant of node j

                    descendents.j <- descendents(nodes = j, dag = dag.G.1)
                    if (!is.character(descendents.j)) {
                      # from score take out parent + descendant
                      sc.tmp <- sc[score.cache$children == j, ,drop=FALSE]
                      sc.tmp <- sc.tmp[sc.tmp[, selected.node] == 1, ,drop=FALSE]

                      # remove descendents of j
                      if (is.null(ncol(sc.tmp))) {
                        if (!is.null(descendents.j)) {
                          if (length(descendents.j) == 1) {
                            sc.tmp <- sc.tmp[sc.tmp[descendents.j] == 0]
                          } else {
                            sc.tmp <- sc.tmp[sum(sc.tmp[descendents.j] == 0) == length(descendents.j)]
                          }
                        }

                      } else {
                        if (!is.null(descendents.j)) {
                          if (length(descendents.j) == 1) {
                            sc.tmp <- sc.tmp[sc.tmp[, descendents.j] == 0, ]
                          } else {
                            sc.tmp <- sc.tmp[rowSums(sc.tmp[, descendents.j,drop=FALSE] == 0) == length(descendents.j), ]
                          }
                        }
                      }

                      dag.G.1[j, ] <- dag.tmp[j, ]

                      if (is.null(ncol(sc.tmp))) {
                        score.G.prime[j] <- logSumExp(sc.tmp[length(sc.tmp)])
                      } else {
                        score.G.prime[j] <- logSumExp(sc.tmp[, ncol(sc.tmp)])
                      }

                    }
                  }
                }

            }
            }
                ############################## Acceptance probability
                n.edges.tilde <- sum(dag.MBR)

                #score.A <-  ( z.i.G.0 / z.star.i.G.prime.0  *  sum(is.finite(score.G)) / sum(is.finite(score.G.prime)))
                #score.A <- n.edges/n.edges.tilde * ((z.star.x.i.M.dot) / (z.star.x.j.M.dot) * (z.x.j.M.cross) / (z.x.i.M.tilde.cross))

                s.proposed <- score.dag(dag.MBR,score.cache,sc)
                s.current <- score.dag(dag.tmp,score.cache,sc)

                A <- min(exp( (s.proposed - s.current) * (n.edges/n.edges.tilde)), 1)

                if(is.nan(A)){A <- 0}

                if (runif(1)<exponent(A, heating)) {
                  rejection <- 0
                  dag.tmp <- dag.MBR

                  score <- score.dag(dag = dag.MBR,bsc.score = score.cache,sc = sc)

                }

    if (is.null(A)) {
        A <- 0
    }
    ############################## Return
    return(list(dag.tmp = dag.tmp, score = score, alpha = A, rejection = rejection))
}  #EOF
