#' GLMM-based Repeatability Using MCMC for Count Data
#' 
#' Calculates repeatability from a generalised linear mixed-effects models fitted by MCMC for count data.
#' 
#' @param y Vector of a response values
#' @param groups Vector of group identities.
#' @param CI Width of the Bayesian credible interval (defaults to 0.95).
#' @param prior List of prior values passed to \link{MCMCglmm} function 
#'        in \pkg{MCMCglmm} (see there for more details). 
#'        Default priors will be used if prior is null.
#' @param verbose Whether or not \link{MCMCglmm} should print MH diagnostics 
#'        are printed to screen. Defaults to FALSE.
#' @param ... Additonal arguements that are passed on to \link{MCMCglmm}
#'        (e.g. length of chain, thinning interval).  
#' 
#' @details Models are fitted using the \link{MCMCglmm} function in
#'          \pkg{MCMCglmm} with \code{poisson} family. Models for binary data 
#'          are fitted with \code{list(R=list(V=1e-10,nu=-1),
#'          G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=25^2)))} unless other
#'          priors are specified in the call.
#' 
#' 
#' @return Returns an object of class rpt that is a a list with the following elements: 
#' \item{datatype}{Type of response (here: count).}
#' \item{method}{Method used to calculate repeatability (here: MCMC).}
#' \item{CI}{Width of the Bayesian credibility interval.}
#' \item{R.link}{Point estimate for repeatability on the link scale, i.e. the 
#'      mode of the posterior distribution.}
#' \item{se.link}{Standard error (\emph{se}) for the repeatability on the link 
#'      scale, i.e. the standard deviation of the posterior distribution. Note
#'      that the distribution might not be symmetrical, in which case \emph{se}
#'      is less informative.}
#' \item{CI.link}{Bayesian credibility interval for the repeatability on the
#'      link scale based on the posterior distribution of \emph{R}.}
#' \item{P.link}{Significance test for the link scale repeatability, returned
#'       as \code{NA}, since the Bayesian approach conflicts with the null 
#'       hypothesis testing.}
#' \item{R.org}{Point estimate for repeatability on the original scale, i.e. 
#'       the mode of the posterior distribution.}
#' \item{se.org}{Standard error (\emph{se}) for repeatability on the original 
#'       scale, i.e. the standard deviation of the posterior distribution. Note 
#'       that the distribution might not be symmetrical, in which case \emph{se}
#'       is less informative.}
#' \item{CI.org}{Bayesian credibility interval for repeatability on the original
#'       scale based on the posterior distribution of \emph{R}.}
#' \item{P.org}{Significance test for the original scale repeatability, 
#'      returned as \code{NA}, since the Bayesian approach conflicts with the 
#'      null hypothesis testing.}
#' \item{R.post}{Named list of MCMC samples form the posterior distributions.
#'       \code{R.link} gives the samples for the link scale repeatability,
#'       \code{R.org} gives the samples for the original scale repeatability.}
#' \item{MCMCpars}{Burnin, length of chain, thinning interval of MCMC chain.}
#' \item{ngroups}{Number of groups.}
#' \item{nobs}{Number of observations.}
#' \item{mod}{Fitted model.}
#'
#'
#' @references 
#' Carrasco, J. L. (2010). \emph{A generalized concordance correlation coefficient based on the variance components generalized linear mixed models with application to overdispersed count data}. Biometrics 66: 897-904.
#'
#' Carrasco, J. L. and Jover, L.  (2005). \emph{Concordance correlation coefficient applied to discrete data}. Statistics in Medicine 24: 4021-4034.
#' 
#' Nakagawa, S. and Schielzeth, H. (2010) \emph{Repeatability for Gaussian and non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
#' 
#' @author Holger Schielzeth  (holger.schielzeth@@ebc.uu.se) & 
#'      Shinichi Nakagawa (shinichi.nakagawa@@otago.ac.nz)
#'      
#' @seealso \link{rpt.poisGLMM.multi}, \link{rpt}, \link{print.rpt}
#' 
#' @examples  
#' # repeatability for female clutch size over two years.
#'    data(BroodParasitism)
#'    attach(BroodParasitism)
#'    (rpt.Host <- rpt.poisGLMM.add(OwnClutches, FemaleID))
#'    detach(BroodParasitism)
#'     
#' # repeatability for male fledgling success
#'    data(Fledglings)
#'    attach(Fledglings)
#'    (rpt.Fledge <- rpt.poisGLMM.add(Fledge, MaleID))
#'    detach(Fledglings) 
#'       
#' @keywords models
#' 
#' @export
#' 
# @importFrom MCMCglmm MCMCglmm
# @importFrom MCMCglmm posterior.mode
# @importFrom coda HPDinterval 

rpt.poisGLMM.add <- function(y, groups, CI=0.95, prior=NULL, verbose=FALSE, ...) {
	# initial checks
	if(any(is.na(y))) warning("missing values in y")
	# preparation
	groups     <- factor(groups)
	# model fitting
	if(is.null(prior)) prior=list(R=list(V=1e-10,nu=-1),G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=25^2)))
	mod        <- MCMCglmm::MCMCglmm(y ~ 1, random=~groups, family="poisson", data=data.frame(y=y,groups=groups), prior=prior, verbose=verbose, ...)
	# ezxtraction of posterior distributions
	var.a      <- mod$VCV[,"groups"]
	var.e      <- mod$VCV[,"units"]
	beta0      <- mod$Sol
	postR.link <- var.a/(var.a + var.e+log(1/exp(beta0)+1))
	EY 		   <- exp(beta0+(var.e+var.a)/2)
	postR.org  <- EY*(exp(var.a)-1)/(EY*(exp(var.e+var.a)-1)+1)
	# point estimate on link and original scale
	R.link     <- MCMCglmm::posterior.mode( postR.link )
	R.org      <- MCMCglmm::posterior.mode( postR.org )
	# credibility interval estimation from paterior distribution
	CI.link    <- coda::HPDinterval(postR.link, CI)[1,]
	CI.org     <- coda::HPDinterval(postR.org, CI)[1,]
	se.link    <- sd(postR.link)
	se.org     <- sd(postR.org)
	# P value equivilents that are not appropriate
	P.link     <- NA
	P.org      <- NA
	# return of results
	res 	   <- list(call=match.call(), datatype="count", method="MCMC", CI=CI,
			      R.link=R.link, se.link=se.link, CI.link=CI.link, P.link=P.link, 
				  R.org = R.org, se.org=se.org, CI.org=CI.org, P.org=P.org,
				  R.post=list(R.link=as.vector(postR.link), R.org=as.vector(postR.org)),
				  MCMCpars = attr(beta0, "mcpar"),
				  ngroups = length(unique(groups)), nobs = length(y),
				  mod = mod)
	class(res) <- "rpt"
	return(res) 
}
