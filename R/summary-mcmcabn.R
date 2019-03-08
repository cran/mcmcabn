############################################################################### summary.mcmcabn.R --- Author : Gilles Kratzer Last modified : 19.02.2019 :

summary.mcmcabn <- function(object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), lag.max = 10, ...) {

    cat("MCMC summary:\n", sep = "")

    cat("Number of Burn in steps: ", (object$burnin), "\n", sep = "")
    cat("Number of MCMC steps: ", (object$iterations), "\n", sep = "")
    cat("Thinning: ", (object$thinning), "\n\n", sep = "")

    cat("Maximum score: ", max(object$scores), "\n", sep = "")
    cat("Empirical mean: ", mean(object$scores), "\n", sep = "")
    cat("Empirical standard deviation: ", sd(object$scores), "\n", sep = "")

    cat("Quantiles of the posterior network score:\n")

    out1 <- matrix(data = quantile(x = object$scores, probs = quantiles), nrow = 1, ncol = length(quantiles), dimnames = list("BN score",
        quantiles))

    print(out1, ...)

    cat("\n\nGlobal acceptance rate: ", 1 - mean(object$rejection), "\n", sep = "")

    out2 <- matrix(data = table(object$method, object$rejection), ncol = 2, dimnames = list(levels(factor(object$method)),
        c("Accepted", "Rejected")))

    print(out2, ...)

    cat("\n\nSample size adjusted for autocorrelation: ", unname(effectiveSize(object$scores)), "\n", sep = "")

    unname(acf(object$scores, lag.max = 10, plot = FALSE))


    cat("\nAutocorrelations by lag:\n")

    out2 <- matrix(data = acf(object$scores, lag.max = 10, plot = FALSE)$acf, nrow = 1, ncol = (lag.max + 1), dimnames = list("acf",
        acf(object$scores, lag.max = 10, plot = FALSE)$lag))

    print(out2, ...)

}  #EOF
