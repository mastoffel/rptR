\name{rpt.remlLMM.adj}
\alias{rpt.remlLMM.adj}
\title{
	LMM-based Repeatability Using REML
}
\description{
	Calculates repeatability from a linear mixed-effects models fitted by REML (restricted maximum likelihood).
}
\usage{
	rpt.remlLMM.adj(formula, grname, data, CI=0.95, nboot=1000, npermut=1000)
}
\arguments{
	\item{formula}{Formula as used e.g. by \link{lmer}. The grouping factor of interest needs to be included as a random effect, e.g. '(1|groups)'. Covariates and additional random effects can be included to estimate adjusted repeatabilities.}
	\item{grname}{A character string or vector of character strings giving the name(s) of the grouping factor(s), for which the repeatability should be estimated. Spelling needs to match the random effect names as given in \code{fromula}.}	
	\item{data}{A dataframe that contains the variables included in the formula argument.}
	\item{CI}{Width of the confidence interval (defaults to 0.95).} 
    \item{nboot}{Number of parametric bootstraps for interval estimation. Defaults to 1000. Larger numbers of permutations give a better asymtotic CI, but may be very time-consuming.}
 	\item{npermut}{Number of permutations used when calculating asymptotic \emph{P} values (defaults to 1000). Currently not in use!}
}
\value{
	Returns an object of class rpt that is a a list with the following elements: 
	\item{datatype}{Response distribution (here: "Gaussian").}
	\item{method}{Method used to calculate repeatability (here: "REML").}
	\item{CI}{Width of the confidence interval.}
	\item{R}{Point estimate for repeatability.}
	\item{se}{Approximate standard error (\emph{se}) for repeatability. Note that the distribution might not be symmetrical, in which case the \emph{se} is less informative.}
	\item{CI.R}{Confidence interval for  repeatability.}
	\item{P}{Approximate \emph{P} value from a significance test based on permutation.}
	\item{R.boot}{Parametric bootstrap samples for \emph{R}.}
	\item{R.permut}{Permutation samples for \emph{R}.}
}
\references{
	Carrasco, J. L. and Jover, L.  (2003). \emph{Estimating the generalized concordance correlation coefficient through variance components}. Biometrics 59: 849-858.
	
	Faraway, J. J. (2006). \emph{Extending the linear model with R}. Boca Raton, FL, Chapman & Hall/CRC.
	
	Nakagawa, S. and Schielzeth, H. (2011) \emph{Repeatability for Gaussian and non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956.
}
\author{
	Holger Schielzeth  (holger.schielzeth@ebc.uu.se) & Shinichi Nakagawa (shinichi.nakagawa@otago.ac.nz)
}
\seealso{
	\link{rpt.mcmcLMM}, \link{print.rpt}, \link{rpt}, \link{rpt.adj}
}
\examples{
\dontrun{
# repeatability estimation for tarsus length - a very high R
data(BodySize)
(rpt.BS <- rpt.remlLMM.adj(Tarsus ~ Sex + (1|BirdID), "BirdID", data=BodySize, nboot=10, npermut=10))
	# reduced number of nboot and npermut iterations

# repeatability estimation for weight (body mass) - a lower R than the previous one
data(BodySize)
(rpt.Weight <- rpt.remlLMM.adj(Weight ~ Sex + (1|BirdID), "BirdID", data=BodySize, nboot=10, npermut=10))
	# reduced number of nboot and npermut iterations
}
}
\keyword{models}