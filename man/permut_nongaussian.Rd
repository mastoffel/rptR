% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internals.R
\name{permut_nongaussian}
\alias{permut_nongaussian}
\title{Permutation function for non-gaussian functions (internal use)}
\usage{
permut_nongaussian(
  permut,
  R_pe,
  formula,
  data,
  dep_var,
  grname,
  npermut,
  parallel,
  ncores,
  link,
  family,
  R,
  rptObj,
  update
)
}
\arguments{
\item{permut}{permutation function which permutes residuals and calculates R}

\item{R_pe}{Function to estimate Repeatabilities and Variances for grouping factors, Residuals, Overdispersion and Fixed effects.}

\item{formula}{lme4 model formula}

\item{data}{data.frame given as original input}

\item{dep_var}{original response variable}

\item{grname}{original grnames vector without Residual or Fixed}

\item{npermut}{number of permutations}

\item{parallel}{boolean}

\item{ncores}{number of cores specified, defaults to NULL}

\item{link}{link function}

\item{family}{respnse family (so far just binomial or poisson)}

\item{R}{point estimate to concetenate with permutations}
}
\description{
Permutation function for non-gaussian functions (internal use)
}
\keyword{internal}
