#' GLMM-based Repeatability Using MCMC for Binomial Data
#' 
#' Calculates repeatability from a generalised linear mixed-effects models 
#' fitted by MCMC for binary and proportion data.
#' 
#' @param data \code{data.frame} containing response and groups variables.
#' @param y If binary, name of response variable in the \code{data.frame}. For proportion data, 
#'        specify two variables: 
#'        The first is the variable for number of successes \emph{m} and the second
#'        variable is \emph{n-m}, where
#'        \emph{m} is the number of successes and \emph{n} is the number of trials.
#' @param groups Name of group variable in the \code{data.frame}.
#' @param CI Width of the Bayesian credible interval (defaults to 0.95)
#' @param prior List of prior values passed to the \link{MCMCglmm} function 
#'        in \pkg{MCMCglmm} (see there for more details). Default priors will be
#'        used if prior is \code{NULL}.
#' @param verbose Whether or not \link{MCMCglmm} should print MH diagnostics 
#'        are printed to screen. Defaults to FALSE.
#' @param ... Additonal arguments that are passed on to \link{MCMCglmm}
#'       (e.g. length of chain, thinning interval).
#'   
#'   
#' @details Models are fitted using the \link{MCMCglmm} function in \pkg{MCMCglmm}.
#'      The categorical family is used for binary data, while the multinomial2 is used for proportion data. 
#'      Models for binary data are fitted with \code{list(R=list(V=1,fix=1),G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=25^2)))} unless other priors are specified in the call.
#'      Models for proportion data are fitted with \code{list(R=list(V=1e-10,nu=-1),G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=25^2)))} unless other priors are specified in the call.
#'   
#' 
#' @return Returns an object of class rpt that is a a list with the following elements: 
#' \item{datatype}{Type of response (here: 'binomial').}
#' \item{method}{Method used to calculate repeatability (here: 'MCMC').}
#' \item{CI}{Width of the Bayesian credibility interval.}
#' \item{R.link}{Point estimate for repeatability on the link scale, i.e. the mode of the posterior distribution.}
#' \item{se.link}{Standard error (\emph{se}) for the repeatability on the link scale, i.e. the standard deviation of the posterior distribution. Note that the distribution might not be symmetrical, in which case \emph{se} is less informative.}
#' \item{CI.link}{Bayesian credibility interval for the intraclass correlation (or repeatability) on the link scale based on the posterior distribution of \emph{R}.}
#' \item{P.link}{Significance test for the link scale repeatability, returned as \code{NA}, since the Bayesian approach conflicts with the null hypothesis testing.}
#' \item{R.org}{Point estimate for repeatability on the original scale, i.e. the mode of the posterior distribution.}
#' \item{se.org}{Standard error (\emph{se}) for repeatability on the original scale, i.e. the standard deviation of the posterior distribution. Note that the distribution might not be symmetrical, in which case \emph{se} is less informative.}
#' \item{CI.org}{Bayesian credibility interval for repeatability on the original scale based on the posterior distribution of \emph{R}.}
#' \item{P.org}{Significance test for the original scale repeatability, returned as \emph{NA}, since the Bayesian approach conflicts with the null hypothesis testing.}
#' \item{R.post}{Named list of MCMC samples form the posterior distributions. \code{R.link} gives the samples for the link scale repeatability, \code{R.org} gives the samples for the original scale repeatability.} 
#' \item{MCMCpars}{Burnin, length of chain, thinning interval of MCMC chain.}
#' \item{ngroups}{Number of groups.}
#' \item{nobs}{Number of observations.}
#' \item{mod}{Fitted model.}
#'
#' @references 
#' Browne, W. J., Subramanian, S. V., et al. (2005). \emph{Variance partitioning in multilevel logistic models that exhibit overdispersion}. Journal of the Royal Statistical Society A 168: 599-613. \cr
#' 
#' Goldstein, H., Browne, W., et al. (2002). \emph{Partitioning variation in multilevel models}. Understanding Statistics 1: 223-231. \cr
#' 
#' Nakagawa, S. and Schielzeth, H. (2010) \emph{Repeatability for Gaussian and non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
#' 
#' @author Holger Schielzeth  (holger.schielzeth@@ebc.uu.se) & 
#'         Shinichi Nakagawa (shinichi.nakagawa@@otago.ac.nz) &
#'         Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'      
#' @seealso \link{rpt.binomGLMM.multi}, \link{rpt}, \link{summary.rpt}, \link{plot.rpt} 
#' 
#' @examples  
#' # repeatability estimations for egg dumping (binary data)
#'      data(BroodParasitism)
#'      EggDump <- subset(BroodParasitism, OwnClutchesBothSeasons == 1, 
#'                        select = c(HostYN, FemaleID))
#'      (rpt.Host <- rpt.binomGLMM.add(data = EggDump, HostYN, FemaleID))
#'      (rpt.BroodPar <- rpt.binomGLMM.add(data = BroodParasitism, cbpYN, FemaleID))
#'      
#' # repeatability estimations for egg dumping (proportion data)
#'      data(BroodParasitism)
#'      ParasitismOR <- subset(BroodParasitism, OwnClutchesBothSeasons == 1, 
#'                             select= c(HostClutches, OwnClutches, FemaleID))
#'      ParasitismOR$parasitised <-  ParasitismOR$OwnClutches -  
#'                                   ParasitismOR$HostClutches 
#'      (rpt.Host <- rpt.binomGLMM.add(data = ParasitismOR, 
#'                                     y = list(HostClutches, parasitised), 
#'                                     groups = FemaleID))
#'  
#' @keywords models
#' 
#' @export
#' 
# @importFrom MCMCglmm MCMCglmm @importFrom MCMCglmm posterior.mode @importFrom coda
# HPDinterval



rpt.binomGLMM.add <- function(data = NULL, y, groups,  CI = 0.95, prior = NULL, verbose = FALSE, ...) {
      
        # data argument should be used
        if (is.null(data)) {
                stop("The data argument needs a data.frame that contains the response (y) and group (groups)")
        }
        rpt.binomGLMM.add_(data, lazyeval::lazy(y), lazyeval::lazy(groups),  CI, prior, verbose , ...)
        
}

#' @export
#' @rdname rpt.binomGLMM.add
#' 
rpt.binomGLMM.add_ <- function(data = NULL, y, groups,  CI = 0.95, prior = NULL, verbose = FALSE, ...) {
    
        if (is.list(lazyeval::lazy_eval(y, data = data))){
                y <- data.frame(do.call(cbind, lazyeval::lazy_eval(y, data = data)))    
        } else {
                y <- lazyeval::lazy_eval(y, data = data)    
        }
        groups <- lazyeval::lazy_eval(groups, data = data)
        
    # data argument check
    if (is.character(y) & ((length(y) == 1) || (length(y) == 2)) & is.character(groups) & 
        (length(groups) == 1)) {
        # y gets vector is one column, stays df if two columns
        ifelse(is.data.frame(data[, y]), y <- as.matrix(data[, y]), y <- data[, y])
        # y <- as.matrix(data[, y])
        groups <- data[, groups]
    }
    # check for equal length y and groups here?
    
    # initial checks
    if (is.null(dim(y))) y <- cbind(y, 1 - y)
    if (any(is.na(y))) warning("missing values in y")
    
    # preparation
    n <- rowSums(y)
    groups <- factor(groups)
    # model fitting # change MCMCglmm with lme4 function plus factor(rownames(data) and
    # family argument --> binomial)
    if (all(n == 1)) {
        if (is.null(prior)) {
            prior <- list(R = list(V = 1, fix = 1), G = list(G1 = list(V = 1, nu = 1, 
                     alpha.mu = 0, alpha.V = 25^2)))
            mod <- MCMCglmm::MCMCglmm(m ~ 1, random = ~groups, data = data.frame(m = y[, 1], 
                   nm = y[, 2], groups = groups), prior = prior, family = "categorical", 
                   verbose = verbose, ...)
        }
    } else {
        if (is.null(prior)) {
            prior <- list(R = list(V = 1e-10, nu = -1), G = list(G1 = list(V = 1, nu = 1, 
                     alpha.mu = 0, alpha.V = 25^2)))
            mod <- MCMCglmm::MCMCglmm(cbind(m, nm) ~ 1, random = ~groups, 
                   data = data.frame(m = y[, 1], nm = y[, 2], groups = groups), 
                   prior = prior, family = "multinomial2", 
                   verbose = verbose, ...)
        }
    }
    # extraction of posterior distributions
    var.a <- mod$VCV[, "groups"]
    var.e <- mod$VCV[, "units"]
    postR.link <- var.a/(var.a + var.e + pi^2/3)
    beta0 <- mod$Sol
    P <- exp(beta0)/(1 + exp(beta0))
    postR.org <- (var.a * P * P/(1 + exp(beta0))^2)/((var.a + var.e) * P * P/(1 + exp(beta0))^2 + 
                 P * (1 - P))
    # point estimate on link and original scale
    R.link <- MCMCglmm::posterior.mode(postR.link)
    R.org <- MCMCglmm::posterior.mode(postR.org)
    # credibility interval estimation from posterior distribution
    CI.link <- coda::HPDinterval(postR.link, CI)[1, ]
    CI.org <- coda::HPDinterval(postR.org, CI)[1, ]
    se.link <- sd(postR.link)
    se.org <- sd(postR.org)
    # P value equivilents that are not appropriate
    P.link <- NA
    P.org <- NA
    # return of results
    res <- list(call = match.call(), datatype = "binomial", method = "MCMC", CI = CI, 
                R.link = R.link, se.link = se.link, CI.link = CI.link, P.link = P.link, R.org = R.org, 
                se.org = se.org, CI.org = CI.org, P.org = P.org, R.post = list(R.link = as.vector(postR.link), 
                R.org = as.vector(postR.org)), MCMCpars = attr(beta0, "mcpar"), ngroups = length(unique(groups)), 
                nobs = length(y), mod = mod)
    class(res) <- "rpt"
    return(res)
} 
