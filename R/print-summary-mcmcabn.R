############################################################################### summary.mcmcabn.R --- Author : Gilles Kratzer Last modified : 19.02.2019 :

print.summary.mcmcabn <- function(x, ...) {

    cat("MCMC summary:\n", sep = "")

    cat("Number of burn-in steps: ", (x$object$burnin), "\n", sep = "")
    cat("Number of MCMC steps: ", (x$object$iterations), "\n", sep = "")
    cat("Thinning: ", (x$object$thinning), "\n\n", sep = "")

    cat("Maximum score: ", max(x$object$scores), "\n", sep = "")
    cat("Empirical mean: ", mean(x$object$scores), "\n", sep = "")
    cat("Empirical standard deviation: ", format(sd(x$object$scores),...), "\n", sep = "")

    cat("Quantiles of the posterior network score:\n")

    print(x$quant, ...)

    cat("\n\nGlobal acceptance rate: ", format( 1 - mean(x$object$rejection), ...), "\n", sep = "")

    print(x$AR, ...)

    cat("\n\nSample size adjusted for autocorrelation: ",
        format( unname(effectiveSize(x$object$scores)), ...), "\n", sep = "")


    cat("\nAutocorrelations by lag:\n")

    print(x$acf, ...)

    invisible(x)
}  #EOF
