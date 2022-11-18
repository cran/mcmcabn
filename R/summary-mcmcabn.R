############################################################################### summary.mcmcabn.R --- Author : Gilles Kratzer Last modified : 19.02.2019 :

summary.mcmcabn <- function(object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), lag.max = 10, ...) {

    out1 <- matrix(data = quantile(x = object$scores, probs = quantiles), nrow = 1,
                   ncol = length(quantiles), dimnames = list("BN score", quantiles))

    out2 <- matrix(data = table(object$method, object$rejection), ncol = 2,
                   dimnames = list(levels(factor(object$method)), c("Accepted", "Rejected")))

    out3 <- matrix(data = acf(object$scores, lag.max = 10, plot = FALSE)$acf, nrow = 1,
                   ncol = (lag.max + 1), dimnames = list("acf", acf(object$scores, lag.max = 10, plot = FALSE)$lag))

    ans <- list( object = object, quant = out1, AR = out2, acf = out3)
    class(ans) <- "summary.mcmcabn"
    return(ans)

}  #EOF
