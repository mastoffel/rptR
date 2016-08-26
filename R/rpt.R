#' (Adjusted) Repeatability Calculation for Gaussian and Non-Gaussian Data
#' 
#' A wrapper function for (adjusted) repeatability estimation from generalized linear mixed-effects 
#' models fitted by restricted maximum likelihood (REML). Calls 
#' specialised functions depending of the choice of datatype and method.
#' 
#' @param formula Formula as used e.g. by \link{lmer}. The grouping factor(s) of
#'        interest needs to be included as a random effect, e.g. '(1|groups)'.
#'        Covariates and additional random effects can be included to estimate adjusted 
#'        repeatabilities.
#' @param grname A character string or vector of character strings giving the
#'        name(s) of the grouping factor(s), for which the repeatability should
#'        be estimated. Spelling needs to match the random effect names as given in \code{formula} 
#'        and terms have to be set in quotation marks. Add "Residual" or "Overdispersion"  to
#'        the character vector to estimate the respective variances. This is most useful
#'        in combination with \code{ratio = FALSE} to estimate the Residual or Overdispersion
#'        variance. With \code{ratio = TRUE} the overdispersion variance reflects the
#'        non-repeatability.
#' @param data A dataframe that contains the variables included in the \code{formula}
#'        and \code{grname} arguments.
#' @param datatype Character string specifying the data type ('Gaussian', 
#'        'Binary', 'Proportion', 'Poisson'). 
#' @param link Character string specifying the link function. Ignored for 
#'        'Gaussian' datatype.
#' @param CI Width of the required confidence interval between 0 and 1 (defaults to 
#'        0.95).
#' @param nboot Number of parametric bootstraps for interval estimation.
#'        (defaults to 1000). Larger numbers of bootstraps give a better
#'        asymtotic CI, but may be very time-consuming (in particular for non-Gaussian data if some variance component 
#'        is low). Bootstrapping can be switch off by setting \code{nboot = 0}.
#' @param npermut Number of permutations used when calculating asymptotic \emph{P} 
#'        values (defaults to 0). Larger numbers of permutations give a better
#'        asymtotic CI, but may be very time-consuming (in particular for non-Gaussian data if some variance component 
#'        is low). Permutaton tests can be switch off by setting \code{npermut = 0}. 
#' @param parallel If TRUE, bootstraps and permutations will be distributed across multiple cores. 
#' @param ncores Specify number of cores to use for parallelization. On default,
#'        all cores but one are used.
#' @param ratio Defaults to TRUE. If FALSE, the variance(s) of the grouping factor(s) of interest
#'        will be used for all further calculations. The resulting point estimate(s), 
#'        uncertainty interval(s) and significance test(s) therefore refer to the estimated variance
#'        itself rather than to the repeatability (i.e. ratio of variances).
#'   
#'   
#' @details For \code{datatype='Gaussian'} calls function \link{rptGaussian},
#'          for \code{datatype='Poisson'} calls function \link{rptPoisson}, 
#'          for \code{datatype='Binary'} calls function \link{rptBinary}, 
#'          for \code{datatype='Proportion'} calls function \link{rptProportion}.
#' 
#' @return Returns an object of class \code{rpt}. See details for specific functions.
#'
#' @references Nakagawa, S. & Schielzeth, H. (2011) \emph{Repeatability for 
#'      Gaussian and non-Gaussian data: a practical guide for biologists}. 
#'      Biological Reviews 85: 935-956.
#'      
#' @author Holger Schielzeth  (holger.schielzeth@@uni-jena.de), 
#'         Shinichi Nakagawa (s.nakagawa@unsw.edu.au),
#'         Martin Stoffel (martin.adam.stoffel@@gmail.com) 
#'         
#' @seealso \link{rptR}
#' 
#' @examples
#' # load data
#' data(BeetlesBody)
#' data(BeetlesMale)
#' data(BeetlesFemale)
#'
#' #  prepare proportion data
#' BeetlesMale$Dark <- BeetlesMale$Colour
#' BeetlesMale$Reddish <- (BeetlesMale$Colour-1)*-1
#' md <- aggregate(cbind(Dark, Reddish) ~ Population + Container, data=BeetlesMale, FUN=sum)
#' 
#' 
#' # Note: all function calls are run without bootstraps and permutation for speed reasons
#' 
#' # gaussian data
#' rpt(BodyL ~ (1|Population), grname="Population", data=BeetlesBody, 
#' nboot=0, npermut=0, datatype = "Gaussian")
#' 
#' # poisson data
#' rpt(Egg ~ Treatment + (1|Container), grname=c("Container"), data = BeetlesFemale,
#' nboot=0, npermut=0, datatype = "Poisson")
#' 
#' # binary data
#' rpt(Colour ~ (1|Population), grname=c("Population"), 
#' data=BeetlesMale, nboot=0, npermut=0, datatype = "Binary")
#' 
#' # proportion data
#' rpt(cbind(Dark, Reddish) ~ (1|Population), grname=c("Population"), data=md,
#' nboot=0, npermut=0, datatype = "Proportion")
#' 
#' 
#' @keywords models
#' 
#' @export
#' 
rpt <- function(formula, grname, data, datatype = c("Gaussian", "Binomial", "Proportion", 
    "count"), link = c("logit", "probit", "log", "sqrt"), CI = 0.95, nboot = 1000, npermut = 0,
    parallel = FALSE, ncores = NULL, ratio = TRUE) {
        
    if (datatype == "Gaussian") {
            return(rptGaussian(formula, grname, data, CI, nboot, npermut, parallel, ncores, ratio))
    }
    if (datatype == "Binary") {
            return(rptBinary(formula, grname, data, link, CI, nboot, npermut, parallel, ncores, ratio))
    }
    if (datatype == "Proportion") {
            return(rptProportion(formula, grname, data, link, CI, nboot, npermut, parallel, ncores, ratio))
    }
    if (datatype == "Poisson") {
            return(rptPoisson(formula, grname, data, link, CI, nboot, npermut, parallel, ncores, ratio))
    }
} 
