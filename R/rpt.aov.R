#' ANOVA-based Repeatability
#' 
#' Repeatability calculations based on Analysis of Variance (ANOVA).
#' 
#' @param y String specifying response variable or vector of measurements.
#' @param groups String specifying groups variable or vector of group identities (will be converted to a factor).
#' @param data Data frame containing respnse and groups variable.
#' @param CI Width of the confidence interval between 0 and 1 (defaults to 0.95).
#' @param npermut Number of permutations used when calculating asymptotic \emph{P} values (defaults to 1000). 
#'   
#' 
#' @return Returns an object of class rpt that is a a list with the following elements: 
#' \item{datatype}{Response distribution (here: 'Gaussian').}
#' \item{method}{Method used to calculate repeatability (here: 'ANOVA').}
#' \item{R}{Point estimate for the repeatability (denoted as \emph{R}).}
#' \item{se}{Asymptotic standard error for repeatability (ICC) based Becker (1982).}
#' \item{CI.R}{Asymptotic confidence interval for repeatability.}
#' \item{P}{Named vector of two \emph{P} values (significance tests): \code{P.aov} is the \emph{P} value for the ANOVA F test, \code{P.permut} is the permutation based \emph{P} value.}
#' \item{R.permut}{Repeatability \emph{R} estimates for each permutation run.}
#' \item{ngroups}{Number of groups.}
#' \item{nobs}{Number of observations.}
#' \item{mod}{Anova model.}
#'
#' @references 
#' Becker, W. A. (1992) \emph{A manual of quantitative genetics}. 5th edn. Academic Enterprises, Pullman, WA. \cr
#' 
#' Lessells, C. M. and Boag, P. T. (1987) \emph{Unrepeatable repeatabilities: a common mistake}. Auk 104: 116-121. \cr
#' 
#' Nakagawa, S. and Schielzeth, H. (2010) \emph{Repeatability for Gaussian and non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
#' 
#' @author Holger Schielzeth  (holger.schielzeth@@ebc.uu.se) & 
#'         Shinichi Nakagawa (shinichi.nakagawa@@otago.ac.nz) &
#'         Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'         
#' @seealso \link{rpt.corr}, \link{rpt.remlLMM}, \link{rpt.mcmcLMM}, \link{rpt}, \link{summary.rpt}, \link{plot.rpt}
#' 
#' @examples  
#' # repeatability estimation for tarsus length - a very high R
#' # reduced number of npermut iterations
#'      data(BodySize)
#'      (rpt.BS <- rpt.aov('Tarsus', 'BirdID', data = BodySize, npermut=10))    
#'      
#' # repeatability estimation for weight (body mass) - a lower R 
#' # than the previous one
#'      data(BodySize)
#'      (rpt.Weight <- rpt.aov('Weight', 'BirdID', data = BodySize, npermut=10))  
#' 
#' @keywords models
#' 
#' @export


rpt.aov <- function(y, groups, data, CI = 0.95, npermut = 1000) {
    
    # check inputs
    if (is.character(y) & (length(y) == 1) & is.character(groups) & (length(groups) == 
        1)) {
        y <- data[[y]]
        groups <- data[[groups]]
    } else if (!(length(y) == length(groups))) {
        stop("y and groups must have the same length")
    }
    
    # initial checks
    if (length(y) != length(groups)) 
        stop("y and group are not of equal length")
    if (npermut < 1) 
        npermut <- 1
    # preparation
    groups <- factor(groups)
    k <- length(levels(groups))
    N <- length(y)
    n0 <- 1/(k - 1) * (N - sum(table(groups)^2)/N)
    # functions: point estimates of R
    R.pe <- function(y, groups, n0) {
        gm <- mean(y)
        MSa <- sum(tapply(y, groups, function(x) (mean(x) - gm)^2 * length(x)))/(k - 
            1)
        MSw <- sum(tapply(y, groups, function(x) sum((x - mean(x))^2)))/(N - k)
        R <- ((MSa - MSw)/n0)/((MSa - MSw)/n0 + MSw)
        return(R)
    }
    # point estimation according to equations 4 and 5
    R <- R.pe(y, groups, n0)
    # confidence interval estimation according to equation 6 and 7
    se <- sqrt((2 * (N - 1) * (1 - R)^2 * (1 + (n0 - 1) * R)^2)/(n0^2 * (N - k) * (k - 
        1)))
    CI.R <- R + c(1, -1) * qt((1 - CI)/2, k - 1) * se
    # significance test from ANOVA
    P.aov <- anova(lm(y ~ groups))[5][1, 1]
    # significance test by permutation
    permut <- function(y, groups, N, n0) {
        sampy <- sample(y, N)
        return(R.pe(sampy, groups, n0))
    }
    if (npermut > 1) {
        R.permut <- c(R, replicate(npermut - 1, permut(y, groups, N, n0), simplify = TRUE))
        P.permut <- sum(R.permut >= R)/npermut
    } else {
        R.permut <- R
        P.permut <- NA
    }
    # return of results
    res <- list(call = match.call(), datatype = "Gaussian", method = "ANOVA", R = R, 
        se = se, CI = CI, CI.R = CI.R, P = c(P.aov = P.aov, P.permut = P.permut), R.permut = R.permut, 
        ngroups = length(unique(groups)), nobs = length(y), mod = anova(lm(y ~ groups)))
    class(res) <- "rpt"
    return(res)
} 
