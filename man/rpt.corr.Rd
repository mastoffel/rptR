\name{rpt.corr}
\alias{rpt.corr}
\title{
	Correlation-based Repeatability
}
\description{
	Calculates repeatability based on inter-class correlations
}
\usage{
	rpt.corr(y, groups, CI=0.95, nboot=1000, npermut=1000)
}
\arguments{
	\item{y}{Vector of measurements. Missing values are not allowed.}
	\item{groups}{Vector of group identitities (will be converted to a factor). Note that each group identity has to appear exactly twice.}
	\item{CI}{Width of the confidence interval between 0 and 1 (defaults to 0.95).}
	\item{nboot}{Number of bootstrapping runs used when calculating an asymptotic confidence interval (defaults to 1000).}
	\item{npermut}{Number of permutations used when calculating asymptotic \emph{P} values (defaults to 1000).}
}
\value{
	Returns an object of class rpt that is a a list with the following elements: 
	\item{datatype}{Response distribution (here: "Gaussian").}
	\item{method}{Method used to calculate repeatability (here: "corr").}
	\item{R}{Point estimate for repeatability \emph{R}.}
	\item{se}{Asymptotic standard error for repeatability based on non-parametric bootstrapping.}
	\item{CI}{Asymptotic confidence interval for repeatability based on non-parametric bootstrapping.}
	\item{P}{Asymptotic \emph{P} value from a significance test for the intraclass correlation based on permutation.}
	\item{R.permut}{Permutation samples for \emph{R}.}
}
\references{
	Sokal, R. R. and F. J. Rohlf (1995). \emph{Biometry: The principles and practice of statistics in biological research}. New York, W.H. Freeman and Company.

	Nakagawa, S. and Schielzeth, H. (2010) \emph{Repeatability for Gaussian and non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
}
\author{
	Holger Schielzeth  (holger.schielzeth@ebc.uu.se) & Shinichi Nakagawa (shinichi.nakagawa@otago.ac.nz)
}
\seealso{
	\link{rpt.aov}, \link{rpt.remlLMM}, \link{rpt.mcmcLMM}, \link{rpt}, \link{print.rpt} 
}
\examples{
# repeatability for male breeding success on a transformed scale
data(Fledglings)
Fledglings$sqrtFledge <- sqrt(Fledglings$Fledge)
attach(Fledglings)
(rpt.Fledge <- rpt.corr(sqrtFledge, MaleID, nboot=10, npermut=10))  # reduced number of iterations
detach(Fledglings)
}
\keyword{models}

