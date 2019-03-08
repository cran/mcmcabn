REV <- function(n.var, dag.tmp, max.parents, sc, score.cache, score, verbose) {

    rejection <- 1
    A <- 0

    # stor number of edges

    n.edges <- sum(dag.tmp)

    # randomly select one

    if (sum(dag.tmp) != 0)
        {

            selected.edge <- which(x = dag.tmp == 1, arr.ind = TRUE)[sample(x = 1:n.edges, size = 1), ]

            # store current parent set

            parent.set <- dag.tmp[selected.edge[1], ]

            # remove parent set
            dag.M.dot <- dag.tmp
            dag.M.dot[selected.edge[1], ] <- 0

            # store descendents j row i col

            descendents.M.dot.i <- descendents(nodes = unname(selected.edge[2]), dag = dag.M.dot)
            descendents.M.dot.j <- descendents(nodes = unname(selected.edge[1]), dag = dag.M.dot)

            # mark all parent sets of i score that do not include j

            sc.tmp <- sc[score.cache$children == selected.edge[2], ]
            sc.tmp <- sc.tmp[sc.tmp[, selected.edge[1]] == 1, ]

            # remove descendents of i
            if (is.null(ncol(sc.tmp))) {
                if (!is.null(descendents.M.dot.i)) {
                  if (length(descendents.M.dot.i) == 1) {
                    sc.tmp <- sc.tmp[sc.tmp[descendents.M.dot.i] == 0]
                  } else {
                    sc.tmp <- sc.tmp[sum(sc.tmp[descendents.M.dot.i] == 0) == length(descendents.M.dot.i)]
                  }
                }
            } else {
                if (!is.null(descendents.M.dot.i)) {
                  if (length(descendents.M.dot.i) == 1) {
                    sc.tmp <- sc.tmp[sc.tmp[, descendents.M.dot.i] == 0, ]
                  } else {
                    sc.tmp <- sc.tmp[rowSums(sc.tmp[, descendents.M.dot.i] == 0) == length(descendents.M.dot.i), ]
                  }
                }
            }



            # sample new parent set

            if (nrow(sc.tmp) > 0 || is.vector(sc.tmp)) {

                if (is.vector(sc.tmp)) {
                  new.parent.i <- sc.tmp
                  new.parent.i <- new.parent.i[-length(new.parent.i)]
                  # store partition function
                  z.star.x.i.M.dot <- sc.tmp[length(sc.tmp)]
                } else {
                  # !!!exp missing
                  new.parent.i <- sc.tmp[sample(x = 1:nrow(sc.tmp), size = 1, prob = range01(sc.tmp[, ncol(sc.tmp)] -
                    sum(sc.tmp[, ncol(sc.tmp)]))), ]
                  new.parent.i <- new.parent.i[-length(new.parent.i)]
                  # store partition function
                  z.star.x.i.M.dot <- (sum(sc.tmp[, ncol(sc.tmp)]))
                }

                # M direct sum
                dag.M.cross <- dag.M.dot
                dag.M.cross[selected.edge[2], ] <- new.parent.i

                # descendant of i
                descendents.M.cross.i <- descendents(nodes = unname(selected.edge[2]), dag = dag.M.cross)
                descendents.M.cross.j <- descendents(nodes = unname(selected.edge[1]), dag = dag.M.cross)

                # score node j not descendant of M
                sc.tmp <- sc[score.cache$children == selected.edge[1], ]

                # remove descendents of i
                if (!is.null(descendents.M.cross.j)) {
                  if (length(descendents.M.cross.j) == 1) {
                    sc.tmp <- sc.tmp[sc.tmp[, descendents.M.cross.j] == 0, ]
                  } else {
                    sc.tmp <- sc.tmp[rowSums(sc.tmp[, descendents.M.cross.j] == 0) == length(descendents.M.cross.j),
                      ]
                  }
                }

                # sample new parent set
                if (nrow(sc.tmp) > 0 || is.vector(sc.tmp)) {

                  if (is.vector(sc.tmp)) {
                    new.parent.j <- sc.tmp
                    # print(new.parent.i)
                    new.parent.j <- new.parent.j[-length(new.parent.j)]
                    # print(new.parent.i) store partition function
                    z.x.j.M.cross <- sc.tmp[length(sc.tmp)]
                  } else {
                    ## !!! exp missing
                    new.parent.j <- sc.tmp[sample(x = 1:nrow(sc.tmp), size = 1, prob = range01(sc.tmp[, ncol(sc.tmp)] -
                      sum(sc.tmp[, ncol(sc.tmp)]))), ]
                    new.parent.j <- new.parent.j[-length(new.parent.j)]
                    # store partition function
                    z.x.j.M.cross <- (sum(sc.tmp[, ncol(sc.tmp)]))
                  }

                  # M tilde

                  dag.M.tilde <- dag.M.cross
                  dag.M.tilde[selected.edge[1], ] <- new.parent.j

                  n.edges.tilde <- sum(dag.M.tilde)

                  ############################## computing acceptance proba

                  # score node j that do not include i
                  sc.tmp <- sc[score.cache$children == selected.edge[1], ]
                  sc.tmp <- sc.tmp[sc.tmp[, selected.edge[2]] == 1, ]
                  # remove descendents of i
                  if (is.vector(sc.tmp)) {
                    if (!is.null(descendents.M.dot.j)) {
                      if (length(descendents.M.dot.j) == 1) {

                        sc.tmp <- sc.tmp[sc.tmp[descendents.M.dot.j] == 0]  #
                      } else {
                        # print(class(sc.tmp)) print(descendents.M.dot.j)
                        sc.tmp <- sc.tmp[sum(sc.tmp[descendents.M.dot.j] == 0) == length(descendents.M.dot.j)]
                      }
                    }
                  } else {
                    if (!is.null(descendents.M.dot.j)) {
                      if (length(descendents.M.dot.j) == 1) {

                        sc.tmp <- sc.tmp[sc.tmp[, descendents.M.dot.j] == 0, ]
                      } else {
                        # print(class(sc.tmp)) print(descendents.M.dot.j)
                        sc.tmp <- sc.tmp[rowSums(sc.tmp[, descendents.M.dot.j] == 0) == length(descendents.M.dot.j),
                          ]
                      }
                    }
                  }

                  if (is.vector(sc.tmp)) {
                    z.star.x.j.M.dot <- sc.tmp[length(sc.tmp)]
                  } else {
                    z.star.x.j.M.dot <- (sum(sc.tmp[, ncol(sc.tmp)]))
                  }


                  dag.M.tilde.cross <- dag.M.dot
                  dag.M.tilde.cross[selected.edge[1], ] <- parent.set

                  # descendant
                  descendents.M.tilde.cross.i <- descendents(nodes = unname(selected.edge[2]), dag = dag.M.tilde.cross)


                  sc.tmp <- sc[score.cache$children == selected.edge[2], ]
                  if (!is.null(descendents.M.tilde.cross.i)) {
                    if (length(descendents.M.tilde.cross.i) == 1) {
                      sc.tmp <- sc.tmp[sc.tmp[, descendents.M.tilde.cross.i] == 0, ]
                    } else {
                      sc.tmp <- sc.tmp[rowSums(sc.tmp[, descendents.M.tilde.cross.i] == 0) == length(descendents.M.tilde.cross.i),
                        ]
                    }
                  }

                  if (is.vector(sc.tmp)) {
                    z.x.i.M.tilde.cross <- sc.tmp[length(sc.tmp)]
                  } else {
                    z.x.i.M.tilde.cross <- (sum(sc.tmp[, ncol(sc.tmp)]))
                  }

                  ############################## Acceptance probability

                  score.A <- n.edges/n.edges.tilde * exp(z.star.x.i.M.dot - z.star.x.j.M.dot + z.x.j.M.cross - z.x.j.M.cross)
                  A <- min(1, score.A)

                  if (rbinom(n = 1, size = 1, prob = A) == 1) {
                    rejection <- 0
                    dag.tmp <- dag.M.tilde
                    score.M.tilde <- 0
                    for (a in 1:n.var) {

                      sc.tmp <- sc[score.cache$children == a, ]
                      score.M.tilde <- sum(min(sc.tmp[(apply(sc.tmp, 1, function(x) identical(unname(x[1:n.var]), unname(dag.tmp[a,
                        ])))), n.var + 1]), score.M.tilde)
                    }
                    score <- score.M.tilde
                  }
                }
            }
        }  #EOFif


    ############################## Return

    return(list(dag.tmp = dag.tmp, score = score, alpha = A, rejection = rejection))
}  #EOF
