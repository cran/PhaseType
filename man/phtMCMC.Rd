\name{phtMCMC}
\alias{phtMCMC}
\title{MCMC for dense Phase-type models}
\description{
Markov-chain Monte Carlo (MCMC) sampler for Bayesian inference on Phase-type models where data consist solely of time to entering the absorbing state (Bladt et al., 2003).

Consider an m+1 state continuous-time Markov chain (CTMC) where the final (m+1 st) state is absorbing.  Then, given data consisting of first passage times to the last state, this function performs Bayesian inference on the rate parameters of the latent continuous-time Markov chain where the generator is assumed to be dense and all rates independent.
}
\usage{
phtMCMC(x, states, beta, nu, zeta, n, mhit=1, resume=NULL, silent=FALSE)
}
\arguments{
  \item{x}{
    a vector of absorption times (or times at which censoring occurs).
}
  \item{states}{
    an integer describing how many states (ie what dimension) the fitted Phase-type's generator should be.
}
  \item{beta}{
    a vector of length m representing the Dirichlet prior on the starting state of the latent continuous-time Markov chain.  Entries should sum to 1.
}
  \item{nu}{
    a list of the Gamma shape hyper-parameters for the prior of each parameter in the continuous-time Markov chain.  These should match the dense generator matrix filled \bold{column-wise}.
}
  \item{zeta}{
    a list of the Gamma reciprocal scale hyper-parameters for the prior of each parameter in the continuous-time Markov chain.  These should match the dense generator matrix, with one zeta value \bold{per row}.
}
  \item{n}{
    the total number of MCMC iterations to compute.
}
  \item{mhit}{
    the number of Metropolis-Hastings iterations to perform when sampling the latent process.
}
  \item{resume}{
    \code{NULL} indicates a new chain is to be sampled.  Otherwise, passing in an object of class \code{phtMCMC} (as returned by a call to this function) picks up a previously run MCMC chain at its end an continues for another \code{n} iterations.
}
  \item{silent}{
    setting to \code{TRUE} suppresses the feedback about current iteration.  Highly recommended for batch files or else output fills with unnecessary iteration updates.
}
}
\details{
Usage of this function effectively involves specification of the number of states for the CTMC generator matrix, the parameter priors in Gamma form and providing absorption time data.

If you want to specify structure on your generator (for example, to model a stochastic process about which some underlying mechanism is known), then consider using \code{\link{phtMCMC2}} instead.  This function is best illustrated by the example below which is \href{http://www.louisaslett.com/PhaseType}{more fully discussed on this webpage}.
}
\value{
  \code{phtMCMC} returns an object of class \code{"phtMCMC"}.

  An object of class \code{"phtMCMC"} is a list containing at least the following components:
  \item{samples}{an object of class \code{"mcmc"} containing the MCMC samples for each parameter in the model.}
  \item{data}{a vector containing the original data used in the inference.}
  \item{vars}{a list of the distinct variable names in the CTMC generator.}
  \item{T}{the naming scheme for the generator which was fitted.}
  \item{beta}{the Dirichlet prior probability mass function on the starting states.}
  \item{nu}{a list of the prior nu hyper-parameters used when phtMCMC2 was called.}
  \item{zeta}{a list of the prior zeta hyper-parameters used when phtMCMC2 was called.}
  \item{iterations}{the number of MCMC iterations contained in the object.}
  \item{MHit}{the number of iterations used if \code{method="MHRS"}}
}
\references{
Bladt, M., Gonzalez, A. and Lauritzen, S. L. (2003), \sQuote{The estimation of Phase-type related functionals using Markov chain Monte Carlo methods}, \emph{Scandinavian Actuarial Journal} 2003(4), 280-300.
}
\author{
Louis J.M. Aslett \email{louis@maths.tcd.ie} (\url{http://www.louisaslett.com/})
}
\note{
Please feel free to email \email{louis@maths.tcd.ie} with any queries or if you encounter errors when running this function.
}


\seealso{
\code{\link{phtMCMC2}}
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
# Gamma prior: shape hyperparameters (one per matrix element, columnwise)
nu <- c(24, 24, 1, 180, 1, 24, 180, 1, 24)
# Gamma prior: reciprocal scale hyperparameters (one per matrix row)
zeta <- c(16, 16, 16)
# Define dimension of model to fit
n <- 3
# Perform 20 MCMC iterations (fix inner Metropolis-Hastings to one iteration
# since starts in stationarity here).  Do more in practice!
res <- phtMCMC(x, n, dirpi, nu, zeta, 6, mhit=1)
print(res)
plot(res)
}
\keyword{ phase-type }
\keyword{ bayesian }
\keyword{ inference }
\keyword{ mcmc }% __ONLY ONE__ keyword per line
