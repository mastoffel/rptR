#' LMM-based repeatability estimated using MCMC sampling
#' 
#' Calculates repeatability from a linear mixed-effects models fitted by MCMC
#' 
#' @param y String specifying response variable or vector of measurements.
#' @param groups String specifying groups variable or vector of group identities (will be converted to a factor).
#' @param data Data frame containing respnse and groups variable.
#' @param CI Width of the Bayesian credible interval (defaults to 0.95)
#' @param prior List of prior values passed to the \link{MCMCglmm} function 
#'        in \pkg{MCMCglmm} (see there for more details). Default priors will 
#'        be used if prior is \code{NULL}.
#' @param verbose Whether or not \link{MCMCglmm} should print MH diagnostics
#'        are printed to screen. Defaults to FALSE.
#' @param ... Additonal arguements that are passed on to \link{MCMCglmm}
#'        (e.g. length of chain, thinning interval).  
#' 
#' @details Models are fitted using the \link{MCMCglmm} function in \pkg{MCMCglmm}. 
#'          Models are fitted with \code{prior=list(R=list(V=1,n=10e-2), 
#'          G=list(G1=list(V=1,n=10e-2)))} unless other priors are specified in the call.
#' 
#' 
#' @return Returns an object of class rpt that is a a list with the following elements: 
#'  \item{datatype}{Response distribution (here: 'Gaussian').}
#'  \item{method}{Method used to calculate repeatability (intra-class correlation, ICC) (here: 'MCMC').}
#'  \item{CI}{Width of the Bayesian credibility interval.}
#'  \item{R}{Point estimate for repeatability (intra-class correlation, ICC), i.e. the mode of the posterior distribution.}
#'  \item{se}{Standard error (\emph{se}) for repeatability (ICC), i.e. the standard deviation of the posterior distribution. Note that the distribution might not be symmetrical, in which case the se is less informative.}
#'  \item{CI.R}{Bayesian credibility interval for the repeatability (ICC) based on the posterior distribution of \emph{R}.}
#'  \item{P}{Significace test, returned as  \code{NA}, since the Bayesian approach conflicts with the null hypothesis testing.}
#'  \item{R.post}{MCMC samples form the posterior distributions of \emph{R}.}
#'  \item{MCMCpars}{Burnin, length of chain, thinning interval of MCMC chain.}
#'  \item{ngroups}{Number of groups.}
#'  \item{nobs}{Number of observations.}
#'  \item{mod}{Fitted model.}
#'
#'
#' @references 
#' Carrasco, J. L. and Jover, L.  (2003). \emph{Estimating the generalized concordance correlation coefficient through variance components}. Biometrics 59: 849-858.
#' 
#' Nakagawa, S. and Schielzeth, H. (2010) \emph{Repeatability for Gaussian and non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
#' 
#' @author Holger Schielzeth  (holger.schielzeth@@ebc.uu.se) & 
#'         Shinichi Nakagawa (shinichi.nakagawa@@otago.ac.nz) &
#'         Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'      
#' @seealso \link{rpt.remlLMM}, \link{rpt.aov}, \link{rpt.corr}, \link{rpt}, \link{summary.rpt}, \link{plot.rpt}
#' 
#' @examples  
#' # repeatability estimation for tarsus length - a very high R
#' data(BodySize)
#' (rpt.BS <- rpt.mcmcLMM(Tarsus, BirdID, data = BodySize))  
#'     
#' # repeatability estimation for weight (body mass) - a lower R than the 
#' # previous one
#' data(BodySize)
#' (rpt.Weight <- rpt.mcmcLMM(Weight, BirdID, data = BodySize))
#'       
#' @keywords models
#' 
#' @export
# @importFrom MCMCglmm MCMCglmm @importFrom MCMCglmm posterior.mode @importFrom coda
# HPDinterval

rpt.mcmcLMM <- function(y, groups, data = NULL, CI = 0.95, prior = NULL, verbose = FALSE, ...) {
    
    # NSE or SE
    if (!is.null(data)) {
                y <- lazyeval::lazy(y)
                groups <- lazyeval::lazy(groups)
                y <- lazyeval::lazy_eval(y$expr, data)
                groups <- lazyeval::lazy_eval(groups$expr, data)
    }
        
    # initial checks
    if (length(y) != length(groups)) 
        stop("y and group are of unequal length")
    # preparation
    groups <- factor(groups)
    if (is.null(prior)) prior <- list(R = list(V = 1, n = 0.1), G = list(G1 = list(V = 1, n = 0.1)))
    # point estimation according to model 8 and equation 9
    mod <- MCMCglmm::MCMCglmm(y ~ 1, random = ~groups, family = "gaussian", 
                              data = data.frame(y = y, groups = groups), prior = prior, 
                              verbose = verbose, ...)
    var.a <- mod$VCV[, "groups"]
    var.e <- mod$VCV[, "units"]
    postR <- var.a/(var.a + var.e)
    # point estimate
    R <- MCMCglmm::posterior.mode(postR)
    # credibility interval estimation from posterior distribution
    CI.R <- coda::HPDinterval(postR, CI)[1, ]
    se <- sd(postR)
    # 'significance test'
    P <- NA
    res <- list(call = match.call(), datatype = "Gaussian", method = "LMM.MCMC", CI = CI, 
        R = R, CI.R = CI.R, se = se, P = P, R.post = postR, MCMCpars = attr(mod$Sol, 
            "mcpar"), ngroups = length(unique(groups)), nobs = length(y), mod = mod)
    class(res) <- "rpt"
    return(res)
} 
