\name{rpt.poisGLMM.multi}
\alias{rpt.poisGLMM.multi}
\title{
	GLMM-based Repeatability Using PQL Estimation for Count Data
}
\description{
	Calculates repeatability from a generalised linear mixed-effects models fitted by PQL (penalized-quasi likelihood)  estimation for count data.
}
\usage{
	rpt.poisGLMM.multi(y, groups, link=c("log", "sqrt"), CI=0.95, nboot=1000, npermut=1000)
}
\arguments{
	\item{y}{Vector of a response values.}
	\item{groups}{Vector of group identitites.}
	\item{link}{Link function. \code{log} and \code{sqrt} are allowed, defaults to \code{log}.}
	\item{CI}{Width of the confidence interval (defaults to 0.95).} 
	\item{nboot}{Number of parametric bootstraps for interval estimation. Defaults to 1000. Larger numbers of permutations give better asymtotic CI, but may be very time-consuming.}
	\item{npermut}{Number of permutations for a significance testing. Defaults to 1000. Larger numbers of permutations give better asymptotic \emph{P} values, but may be very time-consuming.}
}
\value{
	Returns an object of class rpt that is a a list with the following elements: 
	\item{datatype}{Type of response (here: "count").}
	\item{method}{Method used to calculate repeatability (here: "PQL").}
	\item{link}{Link function used (here: "log" or "sqrt").}
	\item{CI}{Width of the confidence interval.}
	\item{R.link}{Point estimate for repeatability (ICC) \emph{R} on the link scale, i.e. the mode of the posterior distribution}
	\item{se.link}{Standard error  (\emph{se}) for repeatability (ICC) on the link scale, i.e. the standard deviation of the posterior distribution. Note that the distribution might not be symmetrical, in which case the se is less informative.}
	\item{CI.link}{Vonfidence interval for repeatability (ICC) on the link scale based on the posterior distribution of  \emph{R}}
	\item{P.link}{Approximate \emph{P} value from a significance test for the link scale repeatability}
	\item{R.org}{Point estimate for repeatability (ICC)  \emph{R} on the original scale, i.e. the mode of the posterior distribution}
	\item{se.org}{Standard error (\emph{se}) for repeatability (ICC) on the original scale, i.e. the standard deviation of the posterior distribution. Note that the distribution might not be symmetrical, in which case \emph{se} is less informative.}
	\item{CI.org}{Confidence interval for repeatability (ICC) on the original scale based on the posterior distribution of  \emph{R}}
	\item{P.org}{Approximate \emph{P} value from a a significance test for the original scale repeatability }
	\item{omega}{Multiplicative overdispersion parameter.}
	\item{R.boot}{Named list of parametric bootstap samples for \emph{R}. \code{R.link} gives the samples for the link scale repeatability, \code{R.org} gives the samples for the original scale repeatability.} 
	\item{R.permut}{Named list of permutation samples for \emph{R}. \code{R.link} gives the samples for the link scale repeatability, \code{R.org} gives the samples for the original scale repeatability.}
}
\details{
	Models are fitted using the \link{glmmPQL} function in \pkg{MASS} with quasipoisson family. 
}
\note{
	Confidence intervals and standard errors are inappropriate at high repeatabilities (\emph{omega} < 1), because parametric bootstrapping allows only \emph{omega} greater than or equal to  1.
}
\references{
	Carrasco, J. L. (2010). \emph{A generalized concordance correlation coefficient based on the variance components generalized linear mixed models with application to overdispersed count data}. Biometrics 66: 897-904.

	Carrasco, J. L. and Jover, L.  (2005). \emph{Concordance correlation coefficient applied to discrete data}. Statistics in Medicine 24: 4021-4034.
	
	Nakagawa, S. and Schielzeth, H. (2010) \emph{Repeatability for Gaussian and non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
}

\author{
	Holger Schielzeth  (holger.schielzeth@ebc.uu.se) & Shinichi Nakagawa (shinichi.nakagawa@otago.ac.nz)
}
\seealso{
	\link{rpt.poisGLMM.add}, \link{rpt}, \link{print.rpt}
}
\examples{
\dontrun{
# repeatability for female clutch size over two years.
data(BroodParasitism)
attach(BroodParasitism)
(rpt.Host <- rpt.poisGLMM.multi(OwnClutches, FemaleID, nboot=10, npermut=10))  
	# reduced number of nboot and npermut iterations
detach(BroodParasitism)

# repeatability for male fledgling success
data(Fledglings)
attach(Fledglings)
(rpt.Fledge <- rpt.poisGLMM.multi(Fledge, MaleID, nboot=10, npermut=10))  
	# reduced number of nboot and npermut iterations
detach(Fledglings)
}
}
\keyword{models}