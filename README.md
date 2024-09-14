# PhaseType R package :package:
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![license](https://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](https://www.gnu.org/licenses/gpl-2.0.html)
[![metacran version](https://www.r-pkg.org/badges/version/PhaseType)](https://cran.r-project.org/package=PhaseType)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/grand-total/PhaseType)](https://cran.r-project.org/package=PhaseType)
[![PhaseType status badge](https://louisaslett.r-universe.dev/badges/PhaseType)](https://louisaslett.r-universe.dev/PhaseType)

This is a package for working with Phase-type (PHT) distributions in the R programming language.
The entire of the MCMC portion of the code has been written in optimised C for higher performance and very low memory use, whilst being easy to call from wrapper R functions.



## Definition of a Phase-type Distribution

Consider a continuous-time Markov chain (CTMC) on a finite discrete state space of size $n+1$, where one of the states is absorbing.
Without loss of generality the generator of the chain can be written in the form:

$$\mathbf{T} = \left( \begin{array}{cc}
  \mathbf{S} & \mathbf{s} \\
  \mathbf{0}^\mathrm{T} & 0
\end{array} \right)$$

where $\mathbf{S}$ is the $n \times n$ matrix of transition rates between non-absorbing states; $\mathbf{s}$ is an $n$ dimensional vector of absorption rates; and $\mathbf{0}$ is an $n$ dimensional vector of zeros.
We take $\boldsymbol{\pi}$ as the initial state distribution: an $n$ dimensional vector of probabilities $\left(\sum_i \pi_i=1\right)$ such that $\pi_i$ is the probability of the chain starting in state $i$.

Then, we define a *Phase-type distribution* to be the distribution of the time to absorption of the CTMC with generator $\mathbf{T}$, or equivalently as the first passage time to state $n+1$.
Thus, a Phase-type distribution is a positively supported univariate distribution having distribution and density functions:

$$\begin{array}{rcl}
  F_X(x) &=& 1 - \boldsymbol{\pi}^\mathrm{T} \exp\\{x \mathbf{S}\\} \mathbf{e}\\
  f_X(x) &=& \boldsymbol{\pi}^\mathrm{T} \exp\\{x \mathbf{S}\\} \mathbf{s}
\end{array}
\qquad \mbox{for } x \in [0,\infty)$$

where $\mathbf{e}$ is an $n$ dimensional vector of $1$'s; $x$ is the time to absorption (or equivalently first-passage time to state $n+1$); and $\exp\\{x \mathbf{S}\\}$ is the matrix exponential.
We denote that a random variable $X$ is Phase-type distributed with parameters $\boldsymbol{\pi}$ and $\mathbf{T}$ by $X \sim \mathrm{PHT}(\boldsymbol{\pi},\mathbf{T})$.

Note that $\displaystyle \sum_{j=1}^n S_{ij} = -s_i \ \forall\\,i$, so often a Phase-type is defined merely by providing $\mathbf{S}$, $\mathbf{T}$ then being implicitly known.



## Contact

Please feel free to:

* submit suggestions and bug-reports at: <https://github.com/louisaslett/PhaseType/issues>
* compose an e-mail to: <louis.aslett@durham.ac.uk>



## Install

You can install the latest release directly from [CRAN](https://CRAN.R-project.org/package=PhaseType).

```r
install.packages("PhaseType")
```



## Install development version (not recommended)

You can get the very latest version from R-universe:

```r
install.packages("PhaseType", repos = c("https://louisaslett.r-universe.dev", "https://cloud.r-project.org"))
```

Installing directly from [GitHub](https://github.com) is not supported by the `install.packages` command, but if you wish to compile from source then you could use the [devtools](https://CRAN.R-project.org/package=devtools) package:

```r
install.packages("remotes")
remotes::install_github("louisaslett/PhaseType")
```

Under releases, the tree/commit from which CRAN releases were made are recorded,
so historic source can be downloaded from there.



## Citation

If you use this software, please use the following citation:

Aslett, L. J. M. (2012), MCMC for Inference on Phase-type and Masked System Lifetime Models, PhD thesis, Trinity College Dublin.

```bibtex
@phdthesis{Aslett2012,
  title={MCMC for Inference on Phase-type and Masked System Lifetime Models},
  author={Aslett, L. J. M.},
  year={2012},
  school={Trinity College Dublin}
}
```

Thank-you :smiley:
