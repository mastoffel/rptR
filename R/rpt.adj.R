#' (Adjusted) Repeatability Calculation for Gaussian and Non-Gaussian Data
#' 
#' A wrapper function for (adjusted) repeatability calculations. Calls 
#' specialised functions depending of the choice of datatype and method.
#' 
#' @param formula Formula as used e.g. by \link{lmer}. The grouping factor of 
#'   interest needs to be included as a random effect, e.g. '(1|groups)'. 
#'   Covariates and additional random effects can be included to estimate 
#'   adjusted repeatabilities.
#' @param grname A character string or vector of character strings giving the 
#'   name(s) of the grouping factor(s), for which the repeatability should be 
#'   estimated. Spelling needs to match the random effect names as given in 
#'   \code{formula}.
#' @param data A dataframe that contains the variables included in the formula 
#'   argument.
#' @param datatype Character string specifying the data type ("Gaussian", 
#'   "binomial", "proportion", "count"). "binomial" and "proportion" are 
#'   interchangable and call the same functions.
#' @param method character string specifying the method of calculation. Defaults
#'   to "REML" for Gaussian data and to "GLMM.multi" for binomial and count 
#'   data.
#' @param link Character string specifying the link function. Ignored for 
#'   "Gaussian" datatype and for the "GLMM.add" method.
#' @param CI Width of the confidence interval between 0 and 1 (defaults to 
#'   0.95).
#' @param nboot Number of bootstrapping runs used when calculating the asymtotic
#'   confidence interval (defaults to 1000). Ignored for the "GLMM.add", "corr" 
#'   and "ANOVA" methods.
#' @param npermut Number of permutations used when calculating asymtotic 
#'   \emph{P} values (defaults to 1000). Ignored for the "GLMM.add" method.
#'   
#'   
#' @details For \code{datatype="Gaussian"} calls function \link{rpt.remlLMM.adj}
#'   or rpt.mcmcLMM.adj (methods "REML" and "MCMC", respecitvely) (Note that
#'   rpt.mcmcLMM.adj is not yet implemented).
#'   
#' 
#' @return Returns an object of class rpt. See details for specific functions.
#' \item{datatype}{Type of repsonse ("Gaussian", "binomial" or "count").}
#' \item{method}{Method used to calculate repeatability ("REML", "MCMC", "ANOVA",
#'       "corr", "GLMM.add" or "GLMM.multi").}
#' \item{link}{Link functions used (GLMMs only).}
#' \item{CI}{Width of the confidence interval or Bayesian credibility interval.}
#' \item{R}{Point estimate for repeatability.}
#' \item{R.link}{Point estimate for repeatability on link scale (GLMM only).}
#' \item{R.org}{Point estimate for repeatability on original scale (GLMM only).}
#' \item{se}{Standard error (\emph{se}) for repeatability. Note that the 
#'      distribution might not be symmetrical, in which case the se is less 
#'      informative.}
#' \item{se.link}{Standard error (\emph{se}) for repeatability on link scale 
#'      (GLMM only).}
#' \item{se.org}{Standard error (\emph{se}) for repeatability on original scale 
#'      (GLMM only).}
#' \item{CI.R}{Confidence interval or Bayesian credibility interval for the 
#'      repeatability.}
#' \item{CI.link}{Confidence interval or Bayesian credibility interval for 
#'      repeatability on link scale (GLMM only).}
#' \item{CI.org}{Confidence interval or Bayesian credibility interval for 
#'      repeatability on original scale (GLMM only).}
#' \item{P}{Significace test, returned as \emph{NA} for the Bayesian approach 
#'      conflicts with the null hypothesis testing.}
#' \item{P.link}{Significace test for repeatability on link scale, returned 
#'      as \emph{NA} for the Bayesian approach conflicts with the null 
#'      hypothesis testing.}
#' \item{P.org}{Significace test for repeatability on original scale, returned 
#'      as \emph{NA} for the Bayesian approach conflicts with the null 
#'      hypothesis testing.}
#' \item{R.post}{MCMC samples form the posterior distributions of \emph{R}.} 
#' \item{R.boot}{Parametric bootstrap samples for \emph{R}.}
#' \item{R.permut}{Permutation samples for \emph{R}.}
#'
#'
#' @references Nakagawa, S. and Schielzeth, H. (2011) \emph{Repeatability for 
#'      Gaussian and non-Gaussian data: a practical guide for biologists}. 
#'      Biological Reviews 85: 935-956.
#'      
#' @author Holger Schielzeth  (holger.schielzeth@@ebc.uu.se) & 
#'      Shinichi Nakagawa (shinichi.nakagawa@@otago.ac.nz)
#' @seealso \link{rpt}
#' 
#' @examples  
#' # for Gaussian data - correlation-based repeatability for male breeding success on a transformed scale
#'      data(Fledglings)
#'      Fledglings$sqrtFledge <- sqrt(Fledglings$Fledge)
#'      rpt.Fledge <- rpt.adj(sqrtFledge ~ Age + (1|MaleID), 
#'                    "MaleID", data=Fledglings, datatype="Gaussian", 
#'                    method="REML", nboot=10, npermut=10)   # reduced number of nboot and npermut iterations
#'      data(BodySize)
#'      (rpt.Weight <- rpt.adj(Weight ~ Sex + (1|BirdID), "BirdID", data=BodySize, datatype="Gaussian", 
#'                     method="MCMC"))
#' 
#' @keywords models
#' 
#' @export
#' 
rpt.adj <- function(formula, grname, data,
			   datatype=c("Gaussian", "binomial", "proportion", "count"),  
			   method=c("corr", "ANOVA", "REML", "MCMC", "GLMM.add", "GLMM.multi"),  
			   link=c("logit", "probit", "log", "sqrt"),
			   CI=0.95, nboot=1000, npermut=1000) {
	if(datatype=="Gaussian") {
		if(length(method)>1) {
			warning("Linear mixed model fitted by REML used by default. Change using argument 'method', if required ('corr', 'ANOVA', 'REML' and 'MCMC' allowed for Gaussian data).")
			method<-"REML" 
		}
		if (method=="REML")  return(rpt.remlLMM.adj(formula, grname, data, CI=CI, nboot=nboot, npermut=npermut))
		if (method=="MCMC")  warning("Not jet implemented") # return(rpt.mcmcLMM.adj(formula, grname, data, CI=CI))
		if (method=="ANOVA") warning("Not jet implemented") # return(rpt.aov(y, groups, CI=CI, npermut=npermut))	
		if (method=="corr")  warning("Not jet implemented") # return(rpt.corr(y, groups, CI=CI, nboot=nboot, npermut=npermut)) 
	}
	if(datatype=="binomial" | datatype=="proportion") {
		if(length(method)>1) {
			warning("Generalised linear mixed model with multiplicative overdispersion fitted by PQL used by default. Change using argument 'method', if required ('GLMM.add' and 'GLMM.multi' allowed for Binomial data).")
			method<-"GLMM.multi" 
		}
		if (method=="GLMM.multi") warning("Not jet implemented") # return(rpt.binomGLMM.multi(y, groups, link, CI=CI, nboot=nboot, npermut=npermut))
		if (method=="GLMM.add") warning("Not jet implemented") # return(rpt.binomGLMM.add(y, groups, CI=CI))
	}
	if(datatype=="count") {
		if(length(method)>1) {
			warning("Generalised linear mixed model with multiplicative overdispersion fitted by PQL used by default. Change using argument 'method', if required ('GLMM.add' and 'GLMM.multi' allowed for count data).")
			method<-"GLMM.multi"
		}
		if(length(link)>1) link="log"
		if (method=="GLMM.multi") warning("Not jet implemented") # return(rpt.poisGLMM.multi(y, groups, link, CI=CI, nboot=nboot, npermut=npermut)) 
		if (method=="GLMM.add")   warning("Not jet implemented") # return(rpt.poisGLMM.add.adj(formula, grname, data, CI=CI)) 
	} 
}