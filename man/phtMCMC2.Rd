\name{phtMCMC2}
\alias{phtMCMC2}
\title{MCMC for structured Phase-type models}
\description{
Markov-chain Monte Carlo (MCMC) sampler for Bayesian inference on Phase-type models where data consist solely of time to entering the absorbing state (Aslett and Wilson, 2011, and Aslett, 2012).

Consider an m+1 state continuous-time Markov chain (CTMC) where the final (m+1 st) state is absorbing.  Then, given data consisting of first passage times to the last state, this function performs Bayesian inference on the rate parameters of the latent continuous-time Markov chain where the generator has some fixed structure.
}
\usage{
phtMCMC2(x, TT, beta, nu, zeta, n, censored=rep(FALSE, length(x)),
         C=matrix(1.0, nrow=dim(TT)[1], ncol=dim(TT)[2]), method="ECS", mhit=1,
         resume=NULL, silent=FALSE)
}
\arguments{
  \item{x}{
    a vector of absorption times (or times at which censoring occurs).
}
  \item{TT}{
    a matrix describing the structure of the Phase-type distribution to be used in the inference.  Each element of the matrix should either be a string for the variable name, or zero to indicate the transition there is prohibited.  See Details and Examples sections below.
}
  \item{beta}{
    a vector of length m representing the Dirichlet prior on the starting state of the latent continuous-time Markov chain.  Entries should sum to 1.
}
  \item{nu}{
    a list of the Gamma shape hyper-parameters for the prior of each parameter in the continuous-time Markov chain.  The \code{tag = value} arguments to the list creation should be such that \code{tag} matches the parameter name given in \code{TT} and \code{value} is the numeric value of the hyper-parameter.
}
  \item{zeta}{
    a list of the Gamma reciprocal scale hyper-parameters for the prior of each parameter in the continuous-time Markov chain.  The \code{tag = value} arguments to the list creation should be such that \code{tag} matches the parameter name given in \code{TT} and \code{value} is the numeric value of the hyper-parameter.
}
  \item{n}{
    the total number of MCMC iterations to compute.
}
  \item{censored}{
    a vector of \code{TRUE/FALSE}s of the same length as \code{x} indicating whether the corresponding observation is censored.
}
  \item{C}{
    a numeric matrix of the same dimensions as \code{TT} which specifies any constant scalar multiples of the parameters in \code{TT}.  Entries must be non-negative.
}
  \item{method}{
    the sampling method to use for the latent stochastic process:

    \code{method="ECS"} is Exact Conditional Sampling as described in Aslett and Wilson (2011);

    \code{method="MHRS"} is Metropolis-Hastings + Rejection Sampling as in the original work of Bladt et al. (2003).
% UNTESTED SO EXCLUDED FROM DOCS; and \code{method="DCS"} is a novel adaptation of the Direct Conditional Sampling technique described by Hobolth (2009).
}
  \item{mhit}{
    the number of Metropolis-Hastings iterations to perform when sampling the latent process.  Ignored unless \code{method="MHRS"}.
}
  \item{resume}{
    \code{NULL} indicates a new chain is to be sampled.  Otherwise, passing in an object of class \code{phtMCMC} (as returned by a call to this function) picks up a previously run MCMC chain at its end an continues for another \code{n} iterations.
}
  \item{silent}{
    setting to \code{TRUE} suppresses the feedback about current iteration.  Highly recommended for batch files or else output fills with unnecessary iteration updates.
}
}
\details{
Usage of this function effectively involves specification of the structure of an absorbing CTMC generator matrix, the parameter priors in Gamma form and providing absorption time data.

The generator matrix \code{TT} is specified as a matrix containing strings (ie text) naming the parameters in each element of the matrix.  This allows constraints to be imposed (see Aslett & Wilson, 2011) and is best illustrated by the example below which is \href{https://www.louisaslett.com/PhaseType/}{more fully discussed on this webpage}.
}
\value{
  \code{phtMCMC2} returns an object of class \code{"phtMCMC"}.

  An object of class \code{"phtMCMC"} is a list containing at least the following components:
  \item{samples}{an object of class \code{"mcmc"} containing the MCMC samples for each parameter in the model.}
  \item{data}{a vector containing the original data used in the inference.}
  \item{vars}{a list of the distinct variable names detected in the CTMC generator.}
  \item{TT}{the original CTMC generator passed in when phtMCMC2 was called.}
  \item{beta}{the Dirichlet prior probability mass function on the starting states.}
  \item{nu}{a list of the prior nu hyper-parameters used when phtMCMC2 was called.}
  \item{zeta}{a list of the prior zeta hyper-parameters used when phtMCMC2 was called.}
  \item{iterations}{the number of MCMC iterations contained in the object.}
  \item{censored}{the original vector of \code{TRUE/FALSE}s of the same length as \code{x} indicating whether the corresponding observations are censored.}
  \item{method}{a string indicating the method used for sampling of the latent process}
  \item{MHit}{the number of iterations used if \code{method="MHRS"}}
}
\references{
Aslett, L. J. M. and Wilson, S. P. (2011), \sQuote{Markov chain Monte Carlo for inference on Phase-type models}, ISI 2011 Proceedings, (\url{https://www.louisaslett.com/Proceedings/ISI_2011.html}).

Aslett, L. J. M. (2012), \sQuote{MCMC for Inference on Phase-type and Masked System Lifetime Models}, PhD Thesis, Trinity College Dublin (\url{https://www.louisaslett.com/PhD_Thesis.html}).

Bladt, M., Gonzalez, A. and Lauritzen, S. L. (2003), \sQuote{The estimation of Phase-type related functionals using Markov chain Monte Carlo methods}, \emph{Scandinavian Actuarial Journal} 2003(4), 280-300.
}
\author{
Louis J.M. Aslett \email{louis.aslett@durham.ac.uk} (\url{https://www.louisaslett.com/})
}
\note{
Please feel free to email \email{louis.aslett@durham.ac.uk} with any queries or if you encounter errors when running this function.
}


\seealso{
\code{\link{phtMCMC}}
}
\examples{
# Some pre-simulated absorption times
x <- c(1.45353415045187, 1.85349532001349, 2.01084961814576, 0.505725921290172,
1.56252630012213, 3.41158665930278, 1.52674487509487, 4.3428662377235,
8.03208018151311, 2.41746547476986, 0.38828086509283, 2.61513815012196,
3.39148865480856, 1.82705817807965, 1.42090953713845, 0.851438991331866,
0.0178808867191894, 0.632198596390046, 0.959910259815998, 1.83344199966323)

# Prior on starting state
dirpi <- c(1, 0, 0)
# Define the structure of the Phase-type generator
TT <- matrix(c(0,"R","R",0,"F",0,0,0,"F",0,0,0,0,"F","F",0), 4)
# Gamma prior: shape hyperparameters (one per model parameter)
nu <- list("R"=180, "F"=24)
# Gamma prior: reciprocal scale hyperparameters (one per model parameter)
zeta <- c("R"=16,"F"=16)
# Perform 20 MCMC iterations.  Do more in practice!!
res <- phtMCMC2(x, TT, dirpi, nu, zeta, 20)
print(res)
\donttest{
plot(res)
}
}
\keyword{ phase-type }
\keyword{ bayesian }
\keyword{ inference }
\keyword{ mcmc }% __ONLY ONE__ keyword per line
