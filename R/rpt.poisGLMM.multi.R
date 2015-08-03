#' GLMM-based Repeatability Using PQL Estimation for Count Data
#' 
#' Calculates repeatability from a generalised linear mixed-effects models 
#' fitted by PQL (penalized-quasi likelihood)  estimation for count data.
#' 
#' @param y Vector of a response values.
#' @param groups Vector of group identities.
#' @param data Data frame containing respnse and groups variable.
#' @param link Link function. \code{log} and \code{sqrt} are allowed, defaults to \code{log}.
#' @param CI Width of the confidence interval (defaults to 0.95).
#' @param nboot Number of parametric bootstraps for interval estimation. Defaults to 1000. 
#'        Larger numbers of permutations give better asymtotic CI, but may be very time-consuming.
#' @param npermut}{Number of permutations for a significance testing. Defaults to 1000.
#'        Larger numbers of permutations give better asymptotic \emph{P} values,
#'        but may be very time-consuming.
#' @param parallel If TRUE, bootstraps will be distributed. 
#' @param ncores Specify number of cores to use for parallelization. On default,
#'        all cores are used. 
#'        
#' @details Models are fitted using the \link{glmmPQL} function in \pkg{MASS} 
#' with quasipoisson family. }\note{ Confidence intervals and standard errors 
#' are inappropriate at high repeatabilities (\emph{omega} < 1), because 
#' parametric bootstrapping allows only \emph{omega} greater than or equal to  1.
#' 
#' 
#' @return Returns an object of class rpt that is a a list with the following elements: 
#' \item{datatype}{Type of response (here: 'count').}
#' \item{method}{Method used to calculate repeatability (here: 'PQL').}
#' \item{link}{Link function used (here: 'log' or 'sqrt').}
#' \item{CI}{Width of the confidence interval.}
#' \item{R.link}{Point estimate for repeatability (ICC) \emph{R} on the link
#'      scale, i.e. the mode of the posterior distribution}
#' \item{se.link}{Standard error  (\emph{se}) for repeatability (ICC) on the 
#'      link scale, i.e. the standard deviation of the posterior distribution.
#'      Note that the distribution might not be symmetrical, in which case the 
#'      se is less informative.}
#' \item{CI.link}{Confidence interval for repeatability (ICC) on the link scale
#'      based on the posterior distribution of  \emph{R}}
#' \item{P.link}{Approximate \emph{P} value from a significance test for the 
#'      link scale repeatability}
#' \item{R.org}{Point estimate for repeatability (ICC)  \emph{R} on the original
#'      scale, i.e. the mode of the posterior distribution}
#' \item{se.org}{Standard error (\emph{se}) for repeatability (ICC) on the 
#'      original scale, i.e. the standard deviation of the posterior distribution.
#'      Note that the distribution might not be symmetrical, in which case
#'      \emph{se} is less informative.}
#' \item{CI.org}{Confidence interval for repeatability (ICC) on the original 
#'      scale based on the posterior distribution of  \emph{R}}
#' \item{P.org}{Approximate \emph{P} value from a a significance test for the 
#'      original scale repeatability }
#' \item{omega}{Multiplicative overdispersion parameter.}
#' \item{R.boot}{Named list of parametric bootstap samples for \emph{R}. 
#'      \code{R.link} gives the samples for the link scale repeatability, 
#'      \code{R.org} gives the samples for the original scale repeatability.} 
#' \item{R.permut}{Named list of permutation samples for \emph{R}. \code{R.link}
#'      gives the samples for the link scale repeatability, \code{R.org} gives 
#'      the samples for the original scale repeatability.}
#' \item{ngroups}{Number of groups.}
#' \item{nobs}{Number of observations.}
#' 
#' @references 
#' Carrasco, J. L. (2010). \emph{A generalized concordance correlation 
#'              coefficient based on the variance components generalized 
#'              linear mixed models with application to overdispersed count data}. 
#'              Biometrics 66: 897-904.
#'
#' Carrasco, J. L. and Jover, L.  (2005). \emph{Concordance correlation coefficient 
#'              applied to discrete data}. Statistics in Medicine 24: 4021-4034.
#' 
#' Nakagawa, S. and Schielzeth, H. (2010) \emph{Repeatability for Gaussian and 
#'              non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
#' 
#' @author Holger Schielzeth  (holger.schielzeth@@ebc.uu.se) & 
#'      Shinichi Nakagawa (shinichi.nakagawa@@otago.ac.nz)
#'      
#' @seealso \link{rpt.poisGLMM.add}, \link{rpt}, \link{print.rpt}
#' 
#' @examples  
#' 
#'        # repeatability for female clutch size over two years.
#'        data(BroodParasitism)
#'        (rpt.Host <- rpt.poisGLMM.multi('OwnClutches', 'FemaleID', 
#'                                         data = BroodParasitism, 
#'                                         nboot=10, npermut=10))  
#'        # reduced number of nboot and npermut iterations
#'        
#'        # repeatability for male fledgling success
#'        data(Fledglings)
#'        (rpt.Fledge <- rpt.poisGLMM.multi('Fledge', 'MaleID', 
#'                                           data = Fledglings, nboot=10, 
#'                                           npermut=10))  
#'        # reduced number of nboot and npermut iterations
#'       
#' @keywords models
#' 
#' @export
#' 
# @importFrom MASS glmmPQL @importFrom lme4 VarCorr
#' 
rpt.poisGLMM.multi <- function(y, groups, data, link = c("log", "sqrt"), CI = 0.95, nboot = 1000, 
    npermut = 1000, parallel = FALSE, ncores = 0) {
    
    if (is.character(y) & ((length(y) == 1) || (length(y) == 2)) & is.character(groups) & 
        (length(groups) == 1)) {
        # y gets vector is one column, stays df if two columns
        ifelse(is.data.frame(data[, y]), y <- as.matrix(data[, y]), y <- data[, y])
        # y <- as.matrix(data[, y])
        groups <- data[, groups]
    }
    
    # initial checks
    if (length(y) != length(groups)) 
        stop("y and group hav to be of equal length")
    if (nboot < 0) 
        nboot <- 0
    if (npermut < 1) 
        npermut <- 1
    if (length(link) > 1) 
        link <- link[1]
    if (link != "log" & link != "sqrt") 
        stop("inappropriate link (has to be 'log' or 'sqrt')")
    if (any(is.na(y))) {
        warning("missing values in y are removed")
        groups <- groups[!is.na(y)]
        y <- y[!is.na(y)]
    }
    # preparation
    groups <- factor(groups)
    N <- length(y)
    k <- length(levels(groups))
    # functions
    pqlglmm.pois.model <- function(y, groups, link, returnR = TRUE) {
        mod <- MASS::glmmPQL(y ~ 1, random = ~1 | groups, family = quasipoisson(link = eval(link)), 
            verbose = FALSE)
        VarComp <- lme4::VarCorr(mod)
        beta0 <- as.numeric(mod$coefficients$fixed)
        omega <- (as.numeric(VarComp[2, 1]))
        var.a <- (as.numeric(VarComp[1, 1]))
        if (link == "log") {
            R.link <- var.a/(var.a + omega * log(1/exp(beta0) + 1))
            EY <- exp(beta0 + var.a/2)
            R.org <- EY * (exp(var.a) - 1)/(EY * (exp(var.a) - 1) + omega)
        }
        if (link == "sqrt") {
            R.link <- var.a/(var.a + omega * 0.25)
            R.org <- NA
        }
        if (returnR) 
            return(list(R.link = R.link, R.org = R.org)) else return(list(beta0 = beta0, omega = omega, var.a = var.a))
    }
    # point estimation according to model 17 equations 18-20
    R <- pqlglmm.pois.model(y, groups, link)
    # confidence interval estimation by parametric bootstrapping nboot is not used in
    # function, but necessary to run parSapply over vector in parallelization
    bootstr <- function(nboot, y, groups, k, N, beta0, var.a, omega, link) {
        groupMeans <- rnorm(k, 0, sqrt(var.a))
        if (link == "log") 
            mu <- exp(beta0 + groupMeans[groups])
        if (link == "sqrt") 
            mu <- (beta0 + groupMeans[groups])^2
        if (omega <= 1) 
            y.boot <- rpois(N, mu) else y.boot <- rnbinom(N, size = (mu/(omega - 1)), mu = mu)
        pqlglmm.pois.model(y.boot, groups, link)
    }
    if (nboot > 0 & parallel == TRUE) {
        mod.ests <- pqlglmm.pois.model(y, groups, link, returnR = FALSE)
        if (ncores == 0) {
            ncores <- parallel::detectCores()
            warning("No core number specified: detectCores() is used to detect the number of \n                        cores on the local machine")
        }
        # start cluster
        cl <- parallel::makeCluster(ncores)
        parallel::clusterExport(cl, "pqlglmm.pois.model", envir = environment())
        # parallel computing
        R.boot <- (parallel::parSapply(cl, 1:nboot, bootstr, y = y, groups = groups, 
            k = k, N = N, beta0 = mod.ests$beta0, var.a = mod.ests$var.a, omega = mod.ests$omega, 
            link = link))
        # transform to list
        R.boot <- lapply(data.frame(t(R.boot)), unlist)
        parallel::stopCluster(cl)
    } else if (nboot > 0 & parallel == FALSE) {
        mod.ests <- pqlglmm.pois.model(y, groups, link, returnR = FALSE)
        R.boot <- replicate(nboot, bootstr(y = y, groups = groups, k = k, N = N, beta0 = mod.ests$beta0, 
            var.a = mod.ests$var.a, omega = mod.ests$omega, link = link), simplify = TRUE)
        R.boot <- list(R.link = as.numeric(unlist(R.boot["R.link", ])), R.org = as.numeric(unlist(R.boot["R.org", 
            ])))
    } else {
        R.boot <- list(R.link = NA, R.org = NA)
    }
    CI.link <- quantile(R.boot$R.link, c((1 - CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
    CI.org <- quantile(R.boot$R.org, c((1 - CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
    se.link <- sd(R.boot$R.link, na.rm = TRUE)
    se.org <- sd(R.boot$R.org, na.rm = TRUE)
    # significance test by randomization
    permut <- function(y, groups, N, link) {
        samp <- sample(1:N, N)
        pqlglmm.pois.model(y, groups[samp], link)
    }
    if (npermut > 1) {
        R.permut <- replicate(npermut - 1, permut(y, groups, N, link), simplify = TRUE)
        R.permut <- list(R.link = c(R$R.link, unlist(R.permut["R.link", ])), R.org = c(R$R.org, 
            unlist(R.permut["R.org", ])))
        P.link <- sum(R.permut$R.link >= R$R.link)/npermut
        P.org <- sum(R.permut$R.org >= R$R.org)/npermut
    } else {
        R.permut <- R
        P.link <- NA
        P.org <- NA
    }
    # return of results
    if (mod.ests$omega < 1) 
        warning("omega < 1, therefore CI limits are unreliable")
    res <- list(call = match.call(), datatype = "count", method = "PQL", link = link, 
        CI = CI, R.link = R$R.link, se.link = se.link, CI.link = CI.link, P.link = P.link, 
        R.org = R$R.org, se.org = se.org, CI.org = CI.org, P.org = P.org, omega = mod.ests$omega, 
        R.boot = list(R.link = R.boot$R.link, R.org = R.boot$R.org), R.permut = list(R.link = R.permut$R.link, 
            R.org = R.permut$R.org), ngroups = length(unique(groups)), nobs = length(y))
    class(res) <- "rpt"
    return(res)
} 
