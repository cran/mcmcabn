% plot.mcmcabn.Rd ---
% Author           : Gilles Kratzer
% Created on :       18.02.2019
% Last modification :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plot.mcmcabn}
\alias{plot.mcmcabn}
\title{Function to plot mcmcabn class objects}

\usage{
\method{plot}{mcmcabn}(x,
     max.score = FALSE,
     \dots)
     }

\arguments{
  \item{x}{object of class mcmcabn.}
  \item{max.score}{logical to plot the cumulative maximum network score.}
  \item{\dots}{arguments to be passed to methods.}
  }

\description{Generic function to plot \code{mcmcabn} objects.
}

\details{The plot function for mcmcabn objects is based on \pkg{ggplot2}, \pkg{ggpubr} and \pkg{cowplot} packages. By default, it returns a trace plot with coloured points when MBR and REV methods have been used. It displays histograms on the right of the densities of (MC)^3, MBR and REV MCMC jumps respectively.}

\author{Gilles Kratzer}

\references{

Plotting ability:
H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

Alboukadel Kassambara (2018). ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.2. https://CRAN.R-project.org/package=ggpubr

Claus O. Wilke (2019). cowplot: Streamlined Plot Theme and Plot Annotations for 'ggplot2'. R package version 0.9.4. https://CRAN.R-project.org/package=cowplot

Data:
Scutari, M. (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1-22. doi:http://dx.doi.org/10.18637/jss.v035.i03.
}


\examples{
## Example from the asia dataset from Lauritzen and Spiegelhalter (1988)
## provided by Scutari (2010)
data("mcmc_run_asia")

# plot the mcmc run
plot(mcmc.out.asia)

# plot cumulative max score
plot(mcmc.out.asia, max.score = TRUE)
}
