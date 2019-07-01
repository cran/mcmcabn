MBR <- function(n.var, dag.tmp, max.parents, sc, score.cache, score, verbose) {

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

    parent.set <- dag.tmp[selected.node, ]

    # remove parent set + co parent
    dag.G.0 <- dag.tmp

    tmp <- which(x = dag.G.0[, selected.node] == 1, arr.ind = TRUE)
    if (!is.integer0(tmp)) {
        dag.G.0[tmp[1], -selected.node] <- 0
    }

    dag.G.0[selected.node, ] <- 0

    # descendant of node i

    descendents.i <- descendents(nodes = selected.node, dag = dag.G.0)


    # from score take out parent + descendant
    sc.tmp <- sc[score.cache$children == selected.node, ]
    if (sum(parent.set) == 1) {
        sc.tmp <- sc.tmp[sc.tmp[, c(as.logical(parent.set), FALSE)] == 0, ]
    }
    if (sum(parent.set) > 1) {
        sc.tmp <- sc.tmp[rowSums(sc.tmp[, c(as.logical(parent.set), FALSE)]) == sum(parent.set), ]
    }

    # remove descendents of i
    if (!is.null(descendents.i)) {
        if (length(descendents.i) == 1) {
            if (is.null(dim(sc.tmp))) {
                sc.tmp <- sc.tmp[sc.tmp[descendents.i] == 0]
            } else {
                sc.tmp <- sc.tmp[sc.tmp[, descendents.i] == 0, ]
            }

        } else {

            if (is.null(dim(sc.tmp))) {
                sc.tmp <- sc.tmp[sum(sc.tmp[descendents.i] == 0) == length(descendents.i)]
            } else {
                sc.tmp <- sc.tmp[rowSums(sc.tmp[, descendents.i] == 0) == length(descendents.i), ]
            }
        }
    }
    # sample new parent set
    if (nrow(sc.tmp) > 0 || is.vector(sc.tmp))
        {

            if (is.vector(sc.tmp)) {
                new.parent.i <- sc.tmp
                new.parent.i <- new.parent.i[-length(new.parent.i)]
                # store partition function
                z.i.G.0 <- sc.tmp[length(sc.tmp)]
            } else {
                # sample a new parent set
                new.parent.i <- sc.tmp[sample(x = 1:nrow(sc.tmp), size = 1, prob = range01(sc.tmp[, ncol(sc.tmp)] -
                  sum(sc.tmp[, ncol(sc.tmp)]))), ]
                new.parent.i <- new.parent.i[-length(new.parent.i)]

                z.i.G.0 <- (sum(sc.tmp[, ncol(sc.tmp)]))
            }
            # dag G1
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
                  if (!is.character(descendents.j)) {


                    # from score take out parent + descendant
                    sc.tmp <- sc[score.cache$children == j, ]
                    sc.tmp <- sc.tmp[sc.tmp[, selected.node] == 1, ]
                    # print(class(sc.tmp)) remove descendents of i

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
                          sc.tmp <- sc.tmp[rowSums(sc.tmp[, descendents.j] == 0) == length(descendents.j), ]
                        }
                      }
                    }


                    if (nrow(sc.tmp) > 0 || is.vector(sc.tmp)) {
                      if (!identical(unname(sc.tmp[length(sc.tmp)]), numeric(0))) {

                        if (is.vector(sc.tmp)) {
                          new.parent.j <- sc.tmp
                          new.parent.j <- new.parent.j[-length(new.parent.j)]

                          # store partition function print(sc.tmp[length(sc.tmp)])

                          score.G[j] <- sc.tmp[length(sc.tmp)]

                        } else {
                          # sample a new parent set
                          new.parent.j <- sc.tmp[sample(x = 1:nrow(sc.tmp), size = 1, prob = range01(sc.tmp[, ncol(sc.tmp)] -
                            sum(sc.tmp[, ncol(sc.tmp)]))), ]
                          new.parent.j <- new.parent.j[-length(new.parent.j)]

                          score.G[j] <- sum(sc.tmp[, ncol(sc.tmp)])
                        }

                        dag.G.1[j, ] <- new.parent.j
                      }
                    }
                  }
                }
            }

            dag.MBR <- dag.G.1

            ############################ complementary inverse move


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
            sc.tmp <- sc[score.cache$children == selected.node, ]
            if (sum(parent.set.G.1) != 0) {
                if (sum(parent.set.G.1) == 1) {
                  sc.tmp <- sc.tmp[sc.tmp[, c(as.logical(parent.set.G.1), FALSE)] == 0, ]
                } else {
                  sc.tmp <- sc.tmp[rowSums(sc.tmp[, c(as.logical(parent.set.G.1), FALSE)]) == sum(parent.set.G.1), ]
                }

                # remove descendents of i
                if (!is.null(descendents.i)) {
                  if (length(descendents.i) == 1) {
                    if (is.null(dim(sc.tmp))) {
                      sc.tmp <- sc.tmp[sc.tmp[descendents.i] == 0]
                    } else {
                      sc.tmp <- sc.tmp[sc.tmp[, descendents.i] == 0, ]
                    }

                  } else {

                    if (is.null(dim(sc.tmp))) {
                      sc.tmp <- sc.tmp[sum(sc.tmp[descendents.i] == 0) == length(descendents.i)]
                    } else {
                      sc.tmp <- sc.tmp[rowSums(sc.tmp[, descendents.i] == 0) == length(descendents.i), ]
                    }
                  }
                }

                # store score
                if (is.null(ncol(sc.tmp))) {
                  z.star.i.G.prime.0 <- sc.tmp[length(sc.tmp)]

                } else {
                  z.star.i.G.prime.0 <- sum(sc.tmp[, ncol(sc.tmp)])
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
                      sc.tmp <- sc[score.cache$children == j, ]
                      sc.tmp <- sc.tmp[sc.tmp[, selected.node] == 1, ]

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
                            sc.tmp <- sc.tmp[rowSums(sc.tmp[, descendents.j] == 0) == length(descendents.j), ]
                          }
                        }
                      }

                      dag.G.1[j, ] <- dag.tmp[j, ]

                      if (is.null(ncol(sc.tmp))) {
                        score.G.prime[j] <- sum(sc.tmp[length(sc.tmp)])
                      } else {
                        score.G.prime[j] <- sum(sc.tmp[, ncol(sc.tmp)])
                      }

                    }
                  }
                }


                ############################## Acceptance probability

                score.MBR <- score.dag(dag = dag.MBR,bsc.score = score.cache,sc = sc)
                #score.A <- exp(z.i.G.0 - z.star.i.G.prime.0 + sum(score.G) - sum(score.G.prime))
                score.A <- exp(- score + score.MBR)


                A <- min(1, score.A)


                if (rbinom(n = 1, size = 1, prob = A) == 1) {
                  rejection <- 0
                  dag.tmp <- dag.MBR
                  #score.MBR <- 0
                  # for (a in 1:n.var) {
                  #   sc.tmp <- sc[score.cache$children == a, ]
                  #   score.MBR <- sum(min(sc.tmp[(apply(sc.tmp, 1, function(x) identical(unname(x[1:n.var]), unname(dag.tmp[a,
                  #     ])))), n.var + 1]), score.MBR)
                  # }
                  score <- score.MBR


                }
            }
        }  #EOIF
    if (is.null(A)) {
        A <- 0
    }
    ############################## Return ##score.MBR <- score
    return(list(dag.tmp = dag.tmp, score = score, alpha = A, rejection = rejection))
}  #EOF
