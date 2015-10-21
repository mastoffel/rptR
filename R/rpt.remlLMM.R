#' LMM-based Repeatability Using REML
#' 
#' Calculates repeatability from a linear mixed-effects models fitted by REML (restricted maximum likelihood).
#' 
#' @param data data.frame containing response and groups variable.
#' @param y Name of response variable in the data.frame. Missing values are not allowed.
#' @param groups Name of group variable in the data.frame.
#' @param CI Width of the confidence interval (defaults to 0.95).
#' @param nboot Number of parametric bootstraps for interval estimation. 
#'        Defaults to 1000. Larger numbers of permutations give a better 
#'        asymtotic CI, but may be very time-consuming.
#' @param npermut Number of permutations used when calculating 
#'        asymptotic \emph{P} values (defaults to 1000).
#' @param parallel If TRUE, bootstraps will be distributed. 
#' @param ncores Specify number of cores to use for parallelization. On default,
#'        all cores are used.
#' 
#' @return 
#' 
#' Returns an object of class rpt that is a a list with the following elements: 
#' \item{call}{Model call.}
#' \item{datatype}{Response distribution (here: 'Gaussian').}
#' \item{method}{Method used to calculate repeatability (here: 'REML').}
#' \item{CI}{Width of the confidence interval.}
#' \item{R}{Point estimate for repeatability.}
#' \item{se}{Approximate standard error (\emph{se}) for repeatability. Note that the distribution might not be symmetrical, in which case the \emph{se} is less informative.}
#' \item{CI.R}{Confidence interval for  repeatability.}
#' \item{P}{Vector of Approximate \emph{P} value from Likelihood-ratio test and Approximate \emph{P} value from a significance test based on permutation.}
#' \item{LRT}{Vector of Likelihood-ratios for the model and the reduced model, and \emph{P} value and degrees of freedom for the Likelihood-ratio test}  
#' \item{R.boot}{Parametric bootstrap samples for \emph{R}.}
#' \item{R.permut}{Permutation samples for \emph{R}.}
#' \item{ngroups}{Number of groups.}
#' \item{nobs}{Number of observations.}
#' \item{mod}{Fitted model.}
#' 
#' @references 
#' Carrasco, J. L. and Jover, L.  (2003). \emph{Estimating the generalized 
#' concordance correlation coefficient through variance components}. Biometrics 59: 849-858.
#'
#' Faraway, J. J. (2006). \emph{Extending the linear model with R}. Boca Raton, FL, Chapman & Hall/CRC.
#' 
#' Nakagawa, S. and Schielzeth, H. (2010) \emph{Repeatability for Gaussian and 
#'              non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
#' 
#' @author Holger Schielzeth  (holger.schielzeth@@ebc.uu.se), 
#'         Shinichi Nakagawa (shinichi.nakagawa@@otago.ac.nz),
#'         Martin A. Stoffel (martin.adam.stoffel@@gmail.com)
#'      
#' @seealso \link{rpt.mcmcLMM}, \link{rpt.aov}, \link{rpt.corr}, \link{print.rpt}, \link{rpt.remlLMM.adj}
#' 
#' @examples  
#' 
#' # repeatability estimation for tarsus length - a very high R
#' data(BodySize)
#' (rpt.BS <- rpt.remlLMM(Tarsus, BirdID, data = BodySize, 
#'                        nboot=10, npermut=10))   
#' # reduced number of nboot and npermut iterations
#'
#' # repeatability estimation for weight (body mass) - a lower R 
#' # than the previous one
#' data(BodySize)
#' (rpt.Weight <- rpt.remlLMM(Weight, BirdID, data = BodySize, 
#'                             nboot=10, npermut=10)) 
#' # reduced number of nboot and npermut iterations
#'       
#'       
#' @keywords models
#' 
#' @export


rpt.remlLMM <- function(data = NULL, y, groups, CI = 0.95, nboot = 1000, npermut = 1000, 
                         parallel = FALSE, ncores = 0) {
        
        # data argument should be used
        if (is.null(data)) {
                stop("The data argument needs a data.frame that contains the response (y) and group (groups)")
        }
        
        rpt.remlLMM_(data = data, lazyeval::lazy(y), lazyeval::lazy(groups), CI = 0.95,
                     nboot, npermut, parallel = FALSE, ncores = 0)
        
}

#' @export 
#' @rdname rpt.remlLMM
rpt.remlLMM_ <- function( data = NULL, y, groups, CI = 0.95, nboot = 1000, npermut = 1000, 
    parallel = FALSE, ncores = 0) {
   
     y <- lazyeval::lazy_eval(y, data = data)
     groups <- lazyeval::lazy_eval(groups, data = data)
     
     # check length equality
     if (!(length(y) == length(groups))) {
                stop("y and groups must have the same length")
      }
    
    # model
    formula <- y ~ 1 + (1 | groups)
    mod <- lme4::lmer(formula)
    
    # checks
    if (nboot < 0) 
        nboot <- 0
    if (npermut < 1) 
        npermut <- 1
    
    # point estimates of R
    R.pe <- function(y, groups) {
        formula <- y ~ 1 + (1 | groups)
        mod.fnc <- lme4::lmer(formula)
        varComps <- lme4::VarCorr(mod.fnc)
        var.a <- as.numeric(varComps)
        var.p <- sum(as.numeric(varComps)) + attr(varComps, "sc")^2
        R <- var.a/var.p
        return(R)
    }
    R <- R.pe(y, groups)
    if (R == 0 & nboot > 0) {
        nboot <- 0
        warning("(One of) the point estimate(s) for the repeatability was exactly zero; parametric bootstrapping has been skipped.")
    }
    
    # confidence interval estimation by parametric bootstrapping
    if (nboot > 0) {
        Ysim <- as.matrix(simulate(mod, nsim = nboot))
    }
    if (nboot > 0 & parallel == TRUE) {
        if (ncores == 0) {
            ncores <- parallel::detectCores()
            warning("No core number specified: detectCores() is used to detect the number of \n cores on the local machine")
        }
        # start cluster
        cl <- parallel::makeCluster(ncores)
        R.boot <- unname(parallel::parApply(cl, Ysim, 2, R.pe, groups = groups))
        parallel::stopCluster(cl)
    } else if (nboot > 0 & parallel == FALSE) {
        R.boot <- unname(apply(Ysim, 2, R.pe, groups = groups))
    } else {
        R.boot <- matrix(rep(NA, length(groups)), nrow = length(groups))
    }
    CI.R <- quantile(R.boot, c((1 - CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
    se <- sd(R.boot)
    
    # significance test by likelihood-ratio-test
    LRT.mod <- logLik(mod)
    LRT.red <- logLik(lm(y ~ 1))
    LRT.D <- as.numeric(-2 * (LRT.red - LRT.mod))
    LRT.df <- 1
    LRT.P <- ifelse(LRT.D <= 0, LRT.df, pchisq(LRT.D, 1, lower.tail = FALSE)/2)
    # significance test by permutation
    permut <- function(formula, groups) {
        groups <- sample(as.character(groups))
        R.pe(y, groups)
    }
    if (npermut > 1) {
        # R.permut <- c(R, replicate(npermut-1, permut(formula, groups), simplify=TRUE))
        R.permut <- c(R, replicate(npermut - 1, permut(formula, groups), simplify = TRUE))
        P.permut <- sum(R.permut >= R)/npermut
    } else {
        R.permut <- R
        P.permut <- NA
    }
    # return of results
    res <- list(call = match.call(), datatype = "Gaussian", method = "LMM.REML", CI = CI, 
        R = R, se = se, CI.R = CI.R, P = c(P.LRT = LRT.P, P.permut = P.permut), LRT = c(LRT.mod = LRT.mod, 
            LRT.red = LRT.red, LRT.D = LRT.D, LRT.df = LRT.df, LRT.P = LRT.P), R.boot = R.boot, 
        R.permut = R.permut, ngroups = length(unique(groups)), nobs = length(y), mod = mod)
    class(res) <- "rpt"
    return(res)
} 
