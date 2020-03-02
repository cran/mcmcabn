% mcmcabn.Rd ---
% Author           : Gilles Kratzer
% Created on :       18.02.2019
% Last modification :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{mcmcabn}
\alias{mcmcabn}
\title{Structural MCMC sampler for DAGs}

\usage{
mcmcabn(score.cache = NULL,
                 score = "mlik",
                 data.dists = NULL,
                 max.parents = 1,
                 mcmc.scheme = c(100,1000,1000),
                 seed = 42,
                 verbose = FALSE,
                 start.dag = NULL,
                 prior.dag = NULL,
                 prior.lambda = NULL,
                 prob.rev = 0.05,
                 prob.mbr = 0.05,
                 heating = 1,
                 prior.choice = 2)
                 }

\arguments{
  \item{score.cache}{output from \link[abn:build_score_cache]{buildscorecache} from the \code{abn} R package.}
  \item{score}{character giving which network score should be used to sample the DAGs landscape.}
\item{data.dists}{a named list giving the distribution for each node in the network, see details.}
  \item{max.parents}{a constant giving the maximum number of parents allowed.}
  \item{mcmc.scheme}{a sampling scheme. It is vector giving in that order: the number of returned DAGS, the number of thinned steps and length of the burn-in phase.}
\item{seed}{a non-negative integer which sets the seed.}
\item{verbose}{extra output, see output for details.}
  \item{start.dag}{a DAG given as a matrix, see details for format, which can be used to provide a starting point for the structural search explicitly. Alternatively, character "random" will select a random DAG as a starting point. Character "hc" will call a hill-climber to select a DAG as a starting point.}
  \item{prior.dag}{user defined prior. It should be given as a matrix where entries range from zero to one. 0.5 is non-informative for the given arc.}
\item{prior.lambda}{hyper parameter representing the strength of belief in the user-defined prior.}
\item{prob.rev}{probability of selecting a new edge reversal.}
\item{prob.mbr}{probability of selecting a Markov blanket resampling move.}
\item{heating}{a real positive number that heats up the chain if between zero and one and apply an exponential decrease scheme if larger than one. The default is one. See details}
  \item{prior.choice}{an integer, 1 or 2, where 1 is a uniform structural prior and 2 uses a weighted prior, see details.}
  }

\description{
This function is a structural Monte Carlo Markov Chain Model Choice (MC)^3 sampler that is equipped with two large scale MCMC moves that are purposed to accelerate chain mixing.
}

\details{The procedure runs a structural Monte Carlo Markov Chain Model Choice (MC)^3 to find the most probable posterior network (DAG). The default algorithm is based on three MCMC moves: edge addition, edge deletion, and edge reversal. This algorithm is known as the (MC)^3. It is known to mix slowly and getting stuck in low probability regions. Indeed, changing of Markov equivalence region often requires multiple MCMC moves. Then large scale MCMC moves are implemented. The user can set the relative frequencies. The new edge reversal move (REV) from Grzegorczyk and Husmeier (2008) and the Markov blanket resampling (MBR) from Su and Borsuk (2016). The classical reversal move depends on the global configuration of the parents and children and fails to propose MCMC jumps that produce valid but very different DAGs in a unique move. The REV move sample globally a new set of parents. The MBR workaround applies the same idea but to the entire Markov blanket of a randomly chosen node.

The classical (MC)^3 is unbiased but inefficient in mixing. The two radical MCMC alternative moves are known to accelerate mixing without introducing biases. Those MCMC moves are computationally expensive. Then low frequencies are advised. The REV move is not necessarily ergotic, then it should not be used alone.

The parameter \code{start.dag} can be: \code{"random"}, \code{"hc"} or user defined. If user select \code{"random"} then a random valid DAG is selected. The routine used favourise low density structure. If \code{"hc"} (for Hill-climber: \link[abn:search_heuristic]{searchHeuristic} then a DAG is selected using 100 different searches with 500 optimization steps. A user defined DAG can be provided. It should be a named square matrix  containing only zeros and ones. The DAG should be valid (i.e. acyclic).

The parameter \code{prior.choice} determines the prior used within each node for a given choice of parent combination. In Koivisto and Sood (2004) p.554, a form of prior is used, which assumes that the prior probability for parent combinations comprising of the same number of parents are all equal. Specifically, that the prior probability for parent set G with cardinality |G| is proportional to 1/[n-1 choose |G|] where there are n total nodes. Note that this favors parent combinations with either very low or very high cardinality, which may not be appropriate. This prior is used when \code{prior.choice=2}. When \code{prior.choice=1} an uninformative prior is used where parent combinations of all cardinalities are equally likely. When \code{prior.choice=3} a user-defined prior is used, defined by \code{prior.dag}. It is given by an adjacency matrix (squared and same size as number of nodes) where entries ranging from zero to one give the user prior belief. An hyperparameter defining the global user belief in the prior is given by \code{prior.lambda}.

MCMC sampler comes with asymptotic statistical guarantees. Therefore it is highly advised to run multiple long enough chains. The burn-in phase length (i.e. throwing away first MCMC iterations) should be adequately chosen.

The argument \code{data.dists} must be a list with named arguments, one for each of the variables in \code{data.df}, where each entry is either \code{"poisson"}, \code{"binomial"}, or \code{"gaussian"}.

The parameter \code{heating} could improve convergence. It should be a real positive number. If smaller than one, it is a tuning parameter which transforms the score by raising it to this power. One is neutral. The smaller, the more probable to accept any move. If larger than one, it indicates the number of returned steps where an exponentially decrease heating scheme is applied. After this number of steps, the \code{heating} parameter is set to one.
}

\value{A list with an entry for the list of sampled DAGs, the list of scores, the acceptance probability, the method used for each MCMC jump, the rejection status for each MCMC jump, the total number of iterations the thinning, the length of burn-in phase, the named list of distribution per node and the heating parameter. The returned object is of class mcmcabn.}

\author{Gilles Kratzer}

\references{

For the implementation of the function:

Kratzer G, Lewis FI, Willi B, Meli ML, Boretti FS, Hofmann-Lehmann R, Torgerson P, Furrer R and Hartnack S (2020) Bayesian Network Modeling Applied to Feline Calicivirus Infection Among Cats in Switzerland. Front. Vet. Sci. 7:73. doi: 10.3389/fvets.2020.00073.

For the new edge reversal:

Grzegorczyk, M., Husmeier, D. (2008). "Improving the structure MCMC sampler for Bayesian networks by introducing a new edge reversal move", Machine Learning, vol. 71(2-3), 265.

For the Markov Blanket resampling move:

Su, C., Borsuk, M. E. (2016). "Improving structure MCMC for Bayesian networks through Markov blanket resampling", The Journal of Machine Learning Research, vol. 17(1), 4042-4061.

For the Koivisto prior:

Koivisto, M. V. (2004). Exact Structure Discovery in Bayesian Networks, Journal of Machine Learning Research, vol 5, 549-573.

For the user defined prior:

Werhli, A. V.,  Husmeier, D. (2007). "Reconstructing gene regulatory networks with Bayesian networks by combining expression data with multiple sources of prior knowledge". Statistical Applications in Genetics and Molecular Biology, 6 (Article 15).

Imoto, S., Higuchi, T., Goto, T., Tashiro, K., Kuhara, S., Miyano, S. (2003). Using Bayesian networks for estimating gene networks from microarrays and biological knowledge. In Proceedings of the European Conference on Computational Biology.

For the asia dataset:

Scutari, M. (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1-22. doi:http://dx.doi.org/10.18637/jss.v035.i03.
}


\examples{
## Example from the asia dataset from Lauritzen and Spiegelhalter (1988)
## provided by Scutari (2010)

# The number of MCMC run is deliberately chosen too small (computing time)
# no thinning (usually not recommended)
# no burn-in (usually not recommended,
# even if not supported by any theoretical arguments)

data("mcmc_run_asia")

# Let us run: 0.03 REV, 0.03 MBR, 0.94 MC3 MCMC jumps
# with a random DAG as starting point

mcmc.out.asia.small <- mcmcabn(score.cache = bsc.compute.asia,
                  score = "mlik",
                  data.dists = dist.asia,
                  max.parents = 2,
                  mcmc.scheme = c(100,0,0),
                  seed = 321,
                  verbose = FALSE,
                  start.dag = "random",
                  prob.rev = 0.03,
                  prob.mbr = 0.03,
                  prior.choice = 2)

summary(mcmc.out.asia.small)

# Soly with MC3 moves

mcmc.out.asia.small <- mcmcabn(score.cache = bsc.compute.asia,
                  score = "mlik",
                  data.dists = dist.asia,
                  max.parents = 2,
                  mcmc.scheme = c(50,0,0),
                  seed = 42,
                  verbose = FALSE,
                  start.dag = "random",
                  prob.rev = 0,
                  prob.mbr = 0,
                  prior.choice = 2)

summary(mcmc.out.asia.small)

# Defining a starting DAG
startDag <- matrix(data = c(0, 0, 0, 1, 0, 0, 0, 0,
                            0, 0, 1, 0, 0, 0, 0, 0,
                            1, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 1, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0),nrow = 8,ncol = 8, byrow = TRUE)

colnames(startDag) <- rownames(startDag) <- names(dist.asia)

# Additionally, let us use the non informative prior
mcmc.out.asia.small <- mcmcabn(score.cache = bsc.compute.asia,
                  score = "mlik",
                  data.dists = dist.asia,
                  max.parents = 2,
                  mcmc.scheme = c(100,0,0),
                  seed = 42,
                  verbose = FALSE,
                  start.dag = startDag,
                  prob.rev = 0,
                  prob.mbr = 0,
                  prior.choice = 1)

summary(mcmc.out.asia.small)

# let us define our very own prior
# we know that there should be a link between Smoking and LungCancer nodes

# uninformative prior matrix
priorDag <- matrix(data = 0.5,nrow = 8,ncol = 8)
# name it
colnames(priorDag) <- rownames(priorDag) <- names(dist.asia)
# parent = smoking; child = LungCancer
priorDag["LungCancer","Smoking"] <- 1

mcmc.out.asia.small <- mcmcabn(score.cache = bsc.compute.asia,
                               score = "mlik",
                               data.dists = dist.asia,
                               max.parents = 2,
                               mcmc.scheme = c(100,0,0),
                               seed = 42,
                               verbose = FALSE,
                               start.dag = startDag,
                               prob.rev = 0,
                               prob.mbr = 0,
                               prior.choice = 3,
                               prior.dag = priorDag)

summary(mcmc.out.asia.small)

# let us improve the convergence rate. The 20 first MCMC moves are performed with an
# heating parameter different than one, afterward, heating is set up to one.

mcmc.out.asia.small <- mcmcabn(score.cache = bsc.compute.asia,
                  score = "mlik",
                  data.dists = dist.asia,
                  max.parents = 2,
                  mcmc.scheme = c(150,0,0),
                  seed = 41242,
                  verbose = FALSE,
                  start.dag = "random",
                  prob.rev = 0.03,
                  prob.mbr = 0.03,
                  prior.choice = 2, heating = 20)

summary(mcmc.out.asia.small)

}
