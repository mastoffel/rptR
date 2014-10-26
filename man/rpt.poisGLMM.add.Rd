\name{rpt.poisGLMM.add}
\alias{rpt.poisGLMM.add}
\title{
	GLMM-based Repeatability Using MCMC for Count Data
}
\description{
	Calculates repeatability from a generalised linear mixed-effects models fitted by MCMC for count data.
}
\usage{
	rpt.poisGLMM.add(y, groups, CI=0.95, prior=NULL, verbose=FALSE, ...)
}
\arguments{
	\item{y}{Vector of a response values.}
	\item{groups}{Vector of group identities.}
	\item{CI}{Width of the Bayesian credible interval (defaults to 0.95).}
	\item{prior}{List of prior values passed to \link{MCMCglmm} function in \pkg{MCMCglmm} (see there for more details). Default priors will be used if prior is null.} 
	\item{verbose}{Whether or not \link{MCMCglmm} should print MH diagnostics are printed to screen. Defaults to FALSE.}
	\item{...}{Additonal arguements that are passed on to \link{MCMCglmm} (e.g. length of chain, thinning interval).}
}
\details{
	Models are fitted using the \link{MCMCglmm} function in \pkg{MCMCglmm} with \code{poisson} family. 
	Models for binary data are fitted with \code{list(R=list(V=1e-10,nu=-1),G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=25^2)))} unless other priors are specified in the call.
}
\value{
	Returns an object of class rpt that is a a list with the following elements: 
	\item{datatype}{Type of response (here: count).}
	\item{method}{Method used to calculate repeatability (here: MCMC).}
	\item{CI}{Width of the Bayesian credibility interval.}
	\item{R.link}{Point estimate for repeatability on the link scale, i.e. the mode of the posterior distribution.}
	\item{se.link}{Standard error (\emph{se}) for the repeatability on the link scale, i.e. the standard deviation of the posterior distribution. Note that the distribution might not be symmetrical, in which case \emph{se} is less informative.}
	\item{CI.link}{Bayesian credibility interval for the repeatability on the link scale based on the posterior distribution of \emph{R}.}
	\item{P.link}{Significance test for the link scale repeatability, returned as \code{NA}, since the Bayesian approach conflicts with the null hypothesis testing.}
	\item{R.org}{Point estimate for repeatability on the original scale, i.e. the mode of the posterior distribution.}
	\item{se.org}{Standard error (\emph{se}) for repeatability on the original scale, i.e. the standard deviation of the posterior distribution. Note that the distribution might not be symmetrical, in which case \emph{se} is less informative.}
	\item{CI.org}{Bayesian credibility interval for repeatability on the original scale based on the posterior distribution of \emph{R}.}
	\item{P.org}{Significance test for the original scale repeatability, returned as \code{NA}, since the Bayesian approach conflicts with the null hypothesis testing.}
	\item{R.post}{Named list of MCMC samples form the posterior distributions. \code{R.link} gives the samples for the link scale repeatability, \code{R.org} gives the samples for the original scale repeatability.} 
}
\references{
	Carrasco, J. L. (2010). \emph{A generalized concordance correlation coefficient based on the variance components generalized linear mixed models with application to overdispersed count data}. Biometrics 66: 897-904.

	Carrasco, J. L. and Jover, L.  (2005). \emph{Concordance correlation coefficient applied to discrete data}. Statistics in Medicine 24: 4021-4034.
	
	Nakagawa, S. and Schielzeth, H. (2010) \emph{Repeatability for Gaussian and non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956.	
}
\author{
	Holger Schielzeth  (holger.schielzeth@ebc.uu.se) & Shinichi Nakagawa (shinichi.nakagawa@otago.ac.nz)
}
\seealso{
	\link{rpt.poisGLMM.multi}, \link{rpt}, \link{print.rpt}
}
\examples{
# repeatability for female clutch size over two years.
data(BroodParasitism)
attach(BroodParasitism)
(rpt.Host <- rpt.poisGLMM.add(OwnClutches, FemaleID))
detach(BroodParasitism)

# repeatability for male fledgling success
data(Fledglings)
attach(Fledglings)
(rpt.Fledge <- rpt.poisGLMM.add(Fledge, MaleID))
detach(Fledglings)
}
\keyword{models}
