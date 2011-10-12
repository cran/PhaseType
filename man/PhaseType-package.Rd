\name{PhaseType-package}
\alias{PhaseType-package}
\alias{PhaseType}
\docType{package}
\title{Phase-type Distributions Toolbox}
\description{
A collection of tools for working with Phase-type Distributions, including sampling methods and both frequentist and Bayesian inference.
}
\details{
\tabular{ll}{
Package: \tab PhaseType\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2011-10-12\cr
License: \tab GPL-2 | GPL-3\cr
LazyLoad: \tab yes\cr
}
%~~ An overview of how to use the package, including the most important ~~
%~~ functions ~~
}
\author{
Louis J. M. Aslett, \email{louis@maths.tcd.ie} (\url{http://www.louisaslett.com})
}
\references{
Aslett, L. J. M. & Wilson, S. P. (2011), Markov chain Monte Carlo for inference on Phase-type models, Technical report, Trinity College Dublin (\emph{pending journal submission}).

Bladt, M., Gonzalez, A. & Lauritzen, S. L. (2003), \sQuote{The estimation of phase-type related functionals using Markov chain Monte Carlo methods}, \emph{Scandinavian Journal of Statistics} \bold{2003}(4), 280-300.
}
\keyword{ phase-type }
\keyword{ pht }
\keyword{ mcmc }
%\seealso{
%%~~ Optional links to other man pages, e.g. ~~
%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
library(actuar)
# Define the S matrix (columnwise)
S <- matrix(c(-3.6, 9.5, 9.5, 1.8, -11.3, 0, 1.8, 0, -11.3), 3)

# Define starting state distribution
pi <- c(1, 0, 0)

# Generate 50 random absorption times from the Phase-type with subgenerator S
# and starting distribution pi, which we will try to infer next
x <- rphtype(50, pi, S)


library(PhaseType)
# FIRST: descriptive model fit (Bladt et al. 2003)
# Prior on starting state
dirpi <- c(1, 0, 0)

# Gamma prior: shape hyperparameters (one per matrix element, columnwise)
nu <- c(24, 24, 1, 180, 1, 24, 180, 1, 24)

# Gamma prior: reciprocal scale hyperparameters (one per matrix row)
zeta <- c(16, 16, 16)

# Define dimension of model to fit
n <- 3

# Perform 20 MCMC iterations (fix inner Metropolis-Hastings to one iteration
# since starts in stationarity here).  Do more in practise!!
res1 <- phtMCMC(x, n, dirpi, nu, zeta, 20, mhit=1)
print(res1)
plot(res1)


# SECOND: mechanistic model fit (Aslett and Wilson 2011)
# Prior on starting state
dirpi <- c(1, 0, 0)

# Define the structure of the Phase-type generator
T <- matrix(c(0,"R","R",0,"F",0,0,0,"F",0,0,0,0,"F","F",0), 4)

# Gamma prior: shape hyperparameters (one per model parameter)
nu <- list("R"=180, "F"=24)

# Gamma prior: reciprocal scale hyperparameters (one per model parameter)
zeta <- c("R"=16,"F"=16)

# Perform 20 MCMC iterations.  Do more in practise!!
res2 <- phtMCMC2(x, T, dirpi, nu, zeta, 20)
print(res2)
plot(res2)
}
