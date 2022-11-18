############################################################################### plot.mcmcabn.R --- Author : Gilles Kratzer Last modified : 19/02/2019 :

plot.mcmcabn <- function(x, max.score = FALSE, ...) {

    # utils::globalVariables(c('X','method' ,'scores')) utils::globalVariables(c('.', '%>%'))
    dta <- data.frame(x[2:4])
    dta$X <- 1:length(x$scores)
    max.score. <- max(x$scores)

    if (!max.score) {
        original_plot <- ggplot(data = dta, aes_string(x = "X", y = "scores")) + geom_line(alpha = 0.5) + geom_hline(yintercept = max.score.,
            linetype = "dashed", color = "red", alpha = 0.8) + geom_text(aes(0, max.score., label = round(max.score.,
            digits = 2), vjust = -0.5), color = "red") + geom_point(data = dta[dta$method %in% c("REV", "MBR"), ], aes_string(color = "factor(method)")) +
            # geom_point(aes(color=as.factor(method)))+
        labs(x = "DAG index", y = "DAG scores", colour = "Methods", title = "Trace plot") + theme_pubr() + scale_color_manual(values = c("#F2C500",
            "#56B4E9"))

        y_density <- axis_canvas(original_plot, axis = "y", coord_flip = TRUE) + geom_density(data = dta, aes_string(x = "scores",
            fill = "factor(method)"), color = NA, alpha = 0.3) + scale_fill_manual(values = c("#F2C500", "#37454B",
            "#56B4E9")) + coord_flip()

        # create the combined plot
        g <- ggdraw(insert_yaxis_grob(plot = original_plot, grob = y_density, position = "right"))
        plot(g)
        invisible(g)
    } else {
        dta$cummax[1] <- dta$scores[1]
        for (i in 2:length(dta$scores)) {
            if (dta$scores[i] > dta$cummax[i - 1]) {
                dta$cummax[i] <- dta$scores[i]
            } else {
                dta$cummax[i] <- dta$cummax[i - 1]
            }
        }
        p <- ggplot(data = dta, aes_string(x = "X", y = "cummax")) + geom_line() + geom_hline(yintercept = max.score.,
            linetype = "dashed", color = "red", alpha = 0.8) + geom_text(aes(0, max.score., label = round(max.score.,
            digits = 2), vjust = -0.5), color = "red") + geom_point(data = dta[dta$method %in% c("REV", "MBR"), ], aes_string(color = "factor(method)")) +
            theme_pubr() + scale_colour_manual(values = c("#F2C500", "#56B4E9")) + labs(x = "DAG index", y = "DAG scores",
            color = "Methods", title = "Trace plot")
        plot(p)
        invisible(p)
    }
}  #EOF
