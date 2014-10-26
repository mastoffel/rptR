\name{rptR}
\alias{rptR}
\title{
	Repeatability Calculation for Gaussian and Non-Gaussian Ddata
}
\description{
	A collection of functions for caluculating point estimates, interval estimates and significance tests of the repeatability (intra-class correlation coefficient) of measurements.
	The function \link{rpt} is a the core functions that calls more specialised functions as required. Specialised functions can also be called directly (see \link{rpt} for details).
	All functions return lists of values. The functions \link{print.rpt} produces summaries in a more appealing format.
}
\note{
	Currently there is are a number of functions that use the two arguments \code{y} and \code{groups} to estimate standard repeatabilities (see e.g. function \link{rpt}). There is another group of functions (e.g. \link{rpt.adj}) that uses the agruments \code{formula} and \code{grname}, which allows estimating adjusted repeatabilities (controlling for fixed effects) and the estimation of multiple variance components simulatneously (multiple random effects). In the long run, the two groups of functions will be merged and will use the more flexible fomula arguments. So far, adjusted repeatabilities are only implemented for Gaussian data using REML estimation.
}
\author{
	Holger Schielzeth  (holger.schielzeth@ebc.uu.se) & Shinichi Nakagawa (shinichi.nakagawa@otago.ac.nz)
}
\references{
	Nakagawa, S. and Schielzeth, H. (2010) \emph{Repeatability for Gaussian and non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
}
\keyword{package}