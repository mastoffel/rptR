\name{rpt.binomGLMM.multi}
\alias{rpt.binomGLMM.multi}
\title{
	GLMM-based Repeatability Using PQL Estimation for Binomial Data
}
\description{
	Calculates repeatability from a generalised linear mixed-effects models fitted by PQL (penalized-quasi likelihood) estimation for binary and proportion data.
}
\usage{
	rpt.binomGLMM.multi(y, groups, link=c("logit", "probit"), CI=0.95, nboot=1000, npermut=1000)
}
\arguments{
	\item{y}{Vector of a response values (for binary data) or a two-column matrix, array or data.frame with colums m, n-m, where m is the number of successes and n the number of trials.}
	\item{groups}{Vector of group identitites.}
	\item{link}{Link function, \code{log} and \code{sqrt} are allowed, defaults to \code{log}.}
	\item{CI}{Width of the confidence interval (defaults to 0.95).} 
	\item{nboot}{Number of parametric bootstraps for interval estimation (defaults to 1000). Larger numbers of permutations give a better asymtotic CI, but may be very time-consuming.}
	\item{npermut}{Number of permutations for significance testing (defaults to 1000). Larger numbers of permutations give better asymtotic \emph{P} values, but may be very time-consuming.}
}
\value{
	Returns an object of class rpt that is a a list with the following elements: 
	\item{datatype}{Type of response (here: "binomial").}
	\item{method}{Method used to calculate repeatability (here: "PQL").}
	\item{link}{Link function used (here: "logit" or "probit").}
	\item{CI}{Width of the confidence interval.}
	\item{R.link}{Point estimate for repeatability on the link scale.}
	\item{se.link}{Standard error  (\emph{se}) for repeatability on the link scale, i.e. the standard deviation of the parametric bootstrap runs. Note that the distribution might not be symmetrical, in which case the \emph{se} is less informative.}
	\item{CI.link}{Confidence interval for repeatability on the link scale based on parametric-boostrapping of \emph{R}.}
	\item{P.link}{Approximate \emph{P} value from a significance test for the link scale repeatability based on randomisation.}
	\item{R.org}{Point estimate for repeatability \emph{R} on the original scale.}
	\item{se.org}{Standard error  (\emph{se}) for repeatability on the original scale, i.e. the standard deviation of the parametric bootstrap runs. Note that the distribution might not be symmetrical, in which case \emph{se} is less informative.}
	\item{CI.org}{Confidence interval for repeatability on the link scale based on parametric-boostrapping of \emph{R}.}
	\item{P.org}{Approximate \emph{P} value from a a significance test for the original scale repeatability based on randomisation. }
	\item{omega}{Multiplicative overdispersion parameter.}
	\item{R.boot}{Named list of parametric bootstap samples for \emph{R}. \code{R.link} gives the samples for the link scale repeatability, \code{R.org} gives the samples for the original scale repeatability.} 
	\item{R.permut}{Named list of permutation samples for \emph{R}. \code{R.link} gives the samples for the link scale repeatability, \code{R.org} gives the samples for the original scale repeatability.}
}
\details{
	Models are fitted using the \link{glmmPQL} function in \pkg{MASS} with the \code{quasibinomial} family (proportion data) or the \code{binomial} family (binary data). 
}
\note{
	Confidence intervals and standard errors are inappropriate at high repeatabilities (\emph{omega} < 1), because parametric bootstrapping allows only \emph{omega} greater than or equal to 1.
}
\references{
	Browne, W. J., Subramanian, S. V., et al. (2005). \emph{Variance partitioning in multilevel logistic models that exhibit overdispersion}. Journal of the Royal Statistical Society A 168: 599-613.
	
	Goldstein, H., Browne, W., et al. (2002). \emph{Partitioning variation in multilevel models} Understanding Statistics 1: 223-231.
	
	Nakagawa, S. and Schielzeth, H. (2011) \emph{Repeatability for Gaussian and non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956.
}
\author{
	Holger Schielzeth  (holger.schielzeth@ebc.uu.se) & Shinichi Nakagawa (shinichi.nakagawa@otago.ac.nz)
}
\seealso{\link{rpt.binomGLMM.add}, \link{rpt}, \link{print.rpt}}
\examples{
\dontrun{
# repeatability estimations for egg dumping (binary data)
data(BroodParasitism)
attach(BroodParasitism)
(rpt.Host <- rpt.binomGLMM.multi(HostYN[OwnClutchesBothSeasons==1], FemaleID[OwnClutchesBothSeasons==1],
	nboot=10, npermut=10)) # low number of nboot and npermut to speed up error checking
(rpt.BroodPar <- rpt.binomGLMM.multi(cbpYN, FemaleID, nboot=10, npermut=10))  
	# low number of nboot and npermut to speed up error checking
detach(BroodParasitism)

# repeatability estimations for egg dumping (proportion data)
data(BroodParasitism)
attach(BroodParasitism)
ParasitisedOR <- cbind(HostClutches, OwnClutches-HostClutches)   
(rpt.Host <- rpt.binomGLMM.multi(ParasitisedOR[OwnClutchesBothSeasons==1,], 
	FemaleID[OwnClutchesBothSeasons==1], nboot=10, npermut=10))  # reduced number of npermut iterations
ParasitismOR <- cbind(cbpEggs, nEggs-cbpEggs)  
zz = which(ParasitismOR[,1]==0 & ParasitismOR[,2]==0) # some rows have entries 0,0 and need to be removed
(rpt.BroodPar <- rpt.binomGLMM.multi(ParasitismOR[-zz,], FemaleID[-zz], nboot=10, npermut=10))   
	# reduced number of npermut iterations
detach(BroodParasitism)
}
}
\keyword{models}