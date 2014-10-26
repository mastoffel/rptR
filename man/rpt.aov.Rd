\name{rpt.aov}
\alias{rpt.aov}
\title{
	ANOVA-based Repeatability
}
\description{
	Repeatability calculations based on Analysis of Variance (ANOVA).
}
\usage{
	rpt.aov(y, groups, CI=0.95, npermut=1000)
}
\arguments{
	\item{y}{Vector of measurements.}
	\item{groups}{Vector of group identities (will be converted to a factor).}
	\item{CI}{Width of the confidence interval between 0 and 1 (defaults to 0.95).}
	\item{npermut}{Number of permutations used when calculating asymptotic \emph{P} values (defaults to 1000).}
}
\value{
	Returns an object of class rpt that is a a list with the following elements: 
	\item{datatype}{Response distribution (here: "Gaussian").}
	\item{method}{Method used to calculate repeatability (here: "ANOVA").}
	\item{R}{Point estimate for the repeatability (denoted as \emph{R}).}
	\item{se}{Asymptotic standard error for repeatability (ICC) based Becker (1982).}
	\item{CI}{Asymptotic confidence interval for repeatability based on non-parametric bootstrapping.}
	\item{P}{Named vector of two \emph{P} values (significance tests): \code{P.aov} is the \emph{P} value for the ANOVA F test, \code{P.permut} is the permutation based \emph{P} value.}
	\item{R.permut}{Repeatability \emph{R} estimates for each permutation run.}
}
\references{
	Becker, W. A. (1992) \emph{A manual of quantitative genetics}. 5th edn. Academic Enterprises, Pullman, WA.
	
	Lessells, C. M. and Boag, P. T. (1987) \emph{Unrepeatable repeatabilities: a common mistake}. Auk 104: 116-121.
	
	Nakagawa, S. and Schielzeth, H. (2010) \emph{Repeatability for Gaussian and non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
}
\author{
	Holger Schielzeth  (holger.schielzeth@ebc.uu.se) & Shinichi Nakagawa (shinichi.nakagawa@otago.ac.nz)
}
\seealso{\link{rpt.corr}, \link{rpt.remlLMM}, \link{rpt.mcmcLMM}, \link{rpt}, \link{print.rpt} }
\examples{
# repeatability estimation for tarsus length - a very high R
data(BodySize)
attach(BodySize)
(rpt.BS <- rpt.aov(Tarsus, BirdID, npermut=10))    # reduced number of npermut iterations
detach(BodySize)

# repeatability estimation for weight (body mass) - a lower R than the previous one
data(BodySize)
attach(BodySize)
(rpt.Weight <- rpt.aov(Weight, BirdID, npermut=10))   # reduced number of npermut iterations
detach(BodySize)
}
\keyword{models}