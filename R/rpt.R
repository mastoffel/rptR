#' Repeatability Estimation for Gaussian and Non-Gaussian Data
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
#'        and terms have to be set in quotation marks. The reseved terms "Residual", 
#'        "Overdispersion" and "Fixed" allow the estimation of oversipersion variance, residual 
#'        variance and variance explained by fixed effects, respectively.
#' @param data A dataframe that contains the variables included in the \code{formula}
#'        and \code{grname} arguments.
#' @param datatype Character string specifying the data type ('Gaussian', 
#'        'Binary', 'Proportion', 'Poisson'). 
#' @param link Character string specifying the link function. Ignored for 
#'        'Gaussian' datatype.
#' @param CI Width of the required confidence interval between 0 and 1 (defaults to 
#'        0.95).
#' @param nboot Number of parametric bootstraps for interval estimation 
#'        (defaults to 1000). Larger numbers of bootstraps give a better
#'        asymtotic CI, but may be time-consuming. Bootstrapping can be switch off by setting 
#'        \code{nboot = 0}. See also \strong{Details} below.
#' @param npermut Number of permutations used when calculating asymptotic p-values 
#'        (defaults to 0). Larger numbers of permutations give a better
#'        asymtotic p-values, but may be time-consuming (in particular when multiple grouping factors
#'        are specified). Permutaton tests can be switch off by setting \code{npermut = 0}. 
#'        See also \strong{Details} below.
#' @param parallel Boolean to express if parallel computing should be applied (defaults to FALSE). 
#'        If TRUE, bootstraps and permutations will be distributed across multiple cores. 
#' @param ncores Specifying the number of cores to use for parallelization. On default,
#'        all but one of the available cores are used.
#' @param ratio Boolean to express if variances or ratios of variance should be estimated. 
#'        If FALSE, the variance(s) are returned without forming ratios. If TRUE (the default) ratios 
#'        of variances (i.e. repeatabilities) are estimated.
#' @param adjusted Boolean to express if adjusted or unadjusted repeatabilities should be estimated. 
#'        If TRUE (the default), the variances explained by fixed effects (if any) will not
#'        be part of the denominator, i.e. repeatabilities are calculated after controlling for 
#'        variation due to covariates. If FALSE, the varianced explained by fixed effects (if any) will
#'        be added to the denominator.
#' @param expect A character string specifying the method for estimating the expectation in Poisson models
#'        with log link and in Binomial models with logit link (in all other cases the agrument is ignored). 
#'        The only valid terms are 'meanobs' and 'latent' (and 'liability for binary and proportion data). 
#'        With the default 'meanobs', the expectation is 
#'        estimated as the mean of the observations in the sample. With 'latent', the expectation is
#'        estimated from estiamtes of the intercept and variances on the link scale. While this is a 
#'        preferred solution, it is susceptible to the distribution of fixed effect covariates and gives 
#'        appropriate results typically only when all covariances are centered to zero. With 'liability' 
#'        estimates follow formulae as presented in Nakagawa & Schielzeth (2010). Liability estimates tend 
#'        to be slightly higher.
#' @param rptObj The output of a rptR function. Can be specified in combination with update = TRUE
#'        to update bootstraps and permutations
#' @param update If TRUE, the rpt object to be updated has to be inputted with the rptObj argument.
#'        The function just updates the permutations and bootstraps, so make sure to specify all other
#'        arguments excactly like for the rpt object specified in rptObj. 
#'   
#' @details 
#' For \code{datatype='Gaussian'} calls function \link{rptGaussian},  
#' for \code{datatype='Poisson'} calls function \link{rptPoisson},   
#' for \code{datatype='Binary'} calls function \link{rptBinary},   
#' for \code{datatype='Proportion'} calls function \link{rptProportion}.  
#'   
#' Confidence intervals and standard errors are estimated by \strong{parametric bootstrapping}. 
#' Under the assumption that the model is specified correctly, the fitted model can be used
#' to generate response values that could potentially be obversed. Differences between the original 
#' data and the simulated response from the fitted model arise from sampling variation. The full model
#' is then fitted to each simuated response vector. The distribution of estimates across all 
#' \code{nboot} replicates represents the design- and model-specific sampling variance and hence 
#' uncertainty of the estimates.
#' 
#' In addition to the likelihood-ratio test, the package uses \strong{permutation tests} for null 
#' hypothesis testing. The general idea is to randomize data under the null hypothesis of no effect 
#' and then test in how many cases the estimates from the model reach or exceed those in the observed 
#' data. In the simplest case, a permutation test randomizes the vector of group identities against 
#' the response vector many times, followed by refitting the model and recalculating the repeatabilities.
#' This provides a null distribution for the case that group identities are unrelated to the response. 
#' However, in more complex models involving multiple random effects and/or fixed effects, such a 
#' procedure will also break the data structure between the grouping factor of interest and other 
#' aspects of the experimental design. Therefore \code{rptR} implements a more robust alternative 
#' which works by fitting a model withouth the grouping factor of interest. It then adds the 
#' randomized residuals to the fitted values of this model, followed by recalculating the repeatability 
#' from the full model. This procedure maintains the general data structure and any effects other 
#' than the grouping effect of interest. The number of permutations can be adjusted with the \code{nperm} argument. 
#' By the logic of a null hypothsis testing, the observed data is one possible (albeit maybe unlikely) 
#' outcome under the null hypothesis. So the observed data is always included as one 'randomization' and 
#' the P value can thus never be lower than \code{1/nperm}, because at least one randomization is as 
#' exteme as the observed data.   
#' 
#' Note also that the \strong{likelihood-ratio test}, since testing variances at the boundary of the 
#' possible parameter range (i.e. against zero), uses a mixture distribution of Chi-square 
#' distrbutions with zero and one degree of freedom as a reference. This ist equivalent to deviding 
#' the P value derived from a Chi-square distribution with one degree of freedom by two.
#' 
#' 
#' @return Returns an object of class \code{rpt}. See specific functions for details.
#'
#' @references Nakagawa, S. & Schielzeth, H. (2010) \emph{Repeatability for 
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
#' BeetlesColour <- aggregate(cbind(Dark, Reddish) ~ Treatment + Population + Container, 
#'      data=BeetlesMale, FUN=sum)
#' 
#' # Note: nboot and npermut are set to 0 for speed reasons. Use larger numbers
#' # for the real analysis.
#' 
#' # gaussian data (example with a single random effect)
#' rpt(BodyL ~ (1|Population), grname="Population", data=BeetlesBody, 
#'      nboot=0, npermut=0, datatype = "Gaussian")
#' 
#' # poisson data (example with two grouping levels and adjusted for fixed effect)
#' rpt(Egg ~ Treatment + (1|Container) + (1|Population), grname=c("Population"), 
#'      data = BeetlesFemale, nboot=0, npermut=0, datatype = "Poisson")
#' 
#' \dontrun{
#' 
#' # binary data (example with estimation of the fixed effect variance)
#' rpt(Colour ~ Treatment + (1|Container) + (1|Population), 
#'      grname=c("Population", "Container", "Fixed"), 
#'      data=BeetlesMale, nboot=0, npermut=0, datatype = "Binary", adjusted = FALSE)
#' 
#' # proportion data (example for the estimation of raw variances, 
#' # including residual and fixed-effect variance)
#' rpt(cbind(Dark, Reddish) ~ Treatment + (1|Population), 
#'      grname=c("Population", "Residual", "Fixed"), data=BeetlesColour,
#'      nboot=0, npermut=0, datatype = "Proportion", ratio=FALSE)
#' 
#' }
#' 
#' @keywords models
#' 
#' @export
#' 
rpt <- function(formula, grname, data, datatype = c("Gaussian", "Binomial", "Proportion", 
    "count"), link = c("logit", "probit", "log", "sqrt"), CI = 0.95, nboot = 1000, npermut = 0,
    parallel = FALSE, ncores = NULL, ratio = TRUE, adjusted = TRUE, expect = "meanobs",
    rptObj = NULL, update = FALSE) {
        
    if (datatype == "Gaussian") {
            out_gaussian <- rptGaussian(formula, grname, data, CI, nboot, npermut, parallel, ncores, ratio, adjusted, rptObj, update)
            out_gaussian$call <- match.call()
            return(out_gaussian)
    }
    if (datatype == "Binary") {
            out_binary <- rptBinary(formula, grname, data, link, CI, nboot, npermut, parallel, ncores, ratio, adjusted, expect, rptObj, update)
            out_binary$call <- match.call()
            return(out_binary)
    }
    if (datatype == "Proportion") {
            out_proportion <- rptProportion(formula, grname, data, link, CI, nboot, npermut, parallel, ncores, ratio, adjusted, expect, rptObj, update)
            out_proportion$call <- match.call()
            return(out_proportion)
    }
    if (datatype == "Poisson") {
            out_poisson <- rptPoisson(formula, grname, data, link, CI, nboot, npermut, parallel, ncores, ratio, adjusted, expect, rptObj, update)
            out_poisson$call <- match.call()
            return(out_poisson)
    }
} 
