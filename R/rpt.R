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
#' @param datatype Character string specifying the data type ('Gaussian', 
#'   'Binary', 'Proportion', 'Poisson'). 
#' @param link Character string specifying the link function. Ignored for 
#'   'Gaussian' datatype.
#' @param CI Width of the confidence interval between 0 and 1 (defaults to 
#'   0.95).
#' @param nboot Number of bootstrapping runs used when calculating the asymtotic
#'   confidence interval (defaults to 1000). 
#' @param npermut Number of permutations used when calculating asymtotic 
#'   \emph{P} values (defaults to 1000). 
#' @param parallel If TRUE, bootstraps will be distributed. 
#' @param ncores Specify number of cores to use for parallelization. On default,
#'        all cores but one are used.
#'   
#' @details For \code{datatype='Gaussian'} calls function \link{rptGaussian},
#'          for \code{datatype='Poisson'} calls function \link{rptPoisson}, etc.
#'   
#' 
#' @return Returns an object of class rpt. See details for specific functions.
#'
#' @references Nakagawa, S. and Schielzeth, H. (2011) \emph{Repeatability for 
#'      Gaussian and non-Gaussian data: a practical guide for biologists}. 
#'      Biological Reviews 85: 935-956.
#'      
#' @author Holger Schielzeth  (holger.schielzeth@@ebc.uu.se),
#'         Shinichi Nakagawa (shinichi.nakagawa@@otago.ac.nz),
#'         Martin A. Stoffel (martin.adam.stoffel@@gmail.com)
#' @seealso \link{rpt}
#' 
#' @examples  
#' 
#' @keywords models
#' 
#' @export
#' 
rpt <- function(formula, grname, data, datatype = c("Gaussian", "Binomial", "Proportion", 
    "count"), link = c("logit", "probit", "log", "sqrt"), CI = 0.95, nboot = 1000, npermut = 1000,
    parallel = FALSE, ncores = NULL) {
        
    if (datatype == "Gaussian") {
            return(rptGaussian(formula, grname, data, CI, nboot, npermut, parallel, ncores))
    }
    if (datatype == "Binary") {
            return(rptBinary(formula, grname, data, link, CI, nboot, npermut, parallel, ncores))
    }
    if (datatype == "Proportion") {
            return(rptProportion(formula, grname, data, link, CI, nboot, npermut, parallel, ncores))
    }
    if (datatype == "Poisson") {
            return(rptPoisson(formula, grname, data, link, CI, nboot, npermut, parallel, ncores))
    }
} 
