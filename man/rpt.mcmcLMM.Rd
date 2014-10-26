\name{rpt.mcmcLMM}
\alias{rpt.mcmcLMM}
\title{
	LMM-based repeatability estimated using MCMC sampling
}
\description{
	Calculates repeatability from a linear mixed-effects models fitted by MCMC
}
\usage{
	rpt.mcmcLMM(y, groups, CI=0.95, prior=NULL, verbose=FALSE, ...)
}
\arguments{
	\item{y}{Vector of a response values}
	\item{groups}{Vector of group identities}
	\item{CI}{Width of the Bayesian credible interval (defaults to 0.95)}
	\item{prior}{List of prior values passed to the \link{MCMCglmm} function in \pkg{MCMCglmm} (see there for more details). Default priors will be used if prior is \code{NULL}.} 
	\item{verbose}{Whether or not \link{MCMCglmm} should print MH diagnostics are printed to screen. Defaults to FALSE.}
	\item{...}{Additonal arguements that are passed on to \link{MCMCglmm} (e.g. length of chain, thinning interval).}
}
\value{
	Returns an object of class rpt that is a a list with the following elements: 
	\item{datatype}{Response distribution (here: "Gaussian").}
	\item{method}{Method used to calculate repeatability (intra-class correlation, ICC) (here: "MCMC").}
	\item{CI}{Width of the Bayesian credibility interval.}
	\item{R}{Point estimate for repeatability (intra-class correlation, ICC), i.e. the mode of the posterior distribution.}
	\item{se}{Standard error (\emph{se}) for repeatability (ICC), i.e. the standard deviation of the posterior distribution. Note that the distribution might not be symmetrical, in which case the se is less informative.}
	\item{CI.R}{Bayesian credibility interval for the repeatability (ICC) based on the posterior distribution of \emph{R}.}
	\item{P}{Significace test, returned as  \code{NA}, since the Bayesian approach conflicts with the null hypothesis testing.}
	\item{R.post}{MCMC samples form the posterior distributions of \emph{R}.} 
}
\details{
	Models are fitted using the \link{MCMCglmm} function in \pkg{MCMCglmm}. 
	Models are fitted with \code{prior=list(R=list(V=1,n=10e-2), G=list(G1=list(V=1,n=10e-2)))} unless other priors are specified in the call.
}
\references{
	Carrasco, J. L. and Jover, L.  (2003). \emph{Estimating the generalized concordance correlation coefficient through variance components}. Biometrics 59: 849-858.
	
	Nakagawa, S. and Schielzeth, H. (2010) \emph{Repeatability for Gaussian and non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
}
\author{
	Holger Schielzeth  (holger.schielzeth@ebc.uu.se) & Shinichi Nakagawa (shinichi.nakagawa@otago.ac.nz)
}
\seealso{
	\link{rpt.remlLMM}, \link{rpt.aov}, \link{rpt.corr}, \link{rpt}, \link{print.rpt}
}
\examples{
# repeatability estimation for tarsus length - a very high R
data(BodySize)
attach(BodySize)
(rpt.BS <- rpt.mcmcLMM(Tarsus, BirdID))  
detach(BodySize)

# repeatability estimation for weight (body mass) - a lower R than the previous one
data(BodySize)
attach(BodySize)
(rpt.Weight <- rpt.mcmcLMM(Weight, BirdID))
detach(BodySize)
}
\keyword{models}