############################################################################### print.mcmcabn.R --- Author : Gilles Kratzer Last modified : 05/02/2019 :

print.mcmcabn <- function(x, ...) {
    cat("Posterior Bayesian network score estimated using MCMC:\n")
    cat("Number of burn-in steps: ", (x$burnin), "\n", sep = "")
    cat("Number of MCMC steps: ", (x$iterations), "\n", sep = "")
    cat("Thinning: ", (x$thinning), "\n\n", sep = "")

    invisible(x)

}  #EOF
