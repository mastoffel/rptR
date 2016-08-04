#' GLMM-based Repeatability Using REML for Gaussian data
#' 
#' Estimates the repeatability from a general linear mixed-effects models fitted by restricted maximum likelihood (REML).
#' @param formula Formula as used e.g. by \link{lmer}. The grouping factor(s) of
#'        interest needs to be included as a random effect, e.g. '(1|groups)'.
#'        Covariates and additional random effects can be included to estimate adjusted 
#'        repeatabilities.
#' @param grname A character string or vector of character strings giving the
#'        name(s) of the grouping factor(s), for which the repeatability should
#'        be estimated. Spelling needs to match the random effect names as given in \code{formula} 
#'        and terms have to be set in quotation marks.
#' @param data A dataframe that contains the variables included in the \code{formula}
#'        and \code{grname} arguments.
#' @param CI Width of the required confidence interval (defaults to 0.95).
#' @param nboot Number of parametric bootstraps for interval estimation.
#'        Defaults to 1000. Larger numbers of bootstraps give a better
#'        asymtotic CI, but may be very time-consuming (in particular of some variance component 
#'        is low). Bootstrapping can be switch off by setting \code{nboot = 0}.
#' @param npermut Number of permutations used when calculating asymptotic \emph{P} 
#'        values (defaults to 1000). Larger numbers of permutations give a better
#'        asymtotic CI, but may be very time-consuming (in particular of some variance component 
#'        is low). Permutaton tests can be switch off by setting \code{npermut = 0}. 
#' @param parallel If TRUE, bootstraps and permutations will be distributed across multiple cores. 
#' @param ncores Specify number of cores to use for parallelization. On default,
#'        all cores but one are used.
#' 
#' @return 
#' Returns an object of class \code{rpt} that is a a list with the following elements: 
#' \item{call}{Function call}
#' \item{datatype}{Response distribution (here: 'Gaussian').}
#' \item{CI}{Coverage of the confidence interval as specified by the \code{CI} argument.}
#' \item{R}{\code{data.frame} with point estimates for repeatabilities for each grouping factor
#'       of interest.}
#' \item{se}{\code{data.frame} with approximate standard errors (\emph{se}) for repeatabilities. 
#'      Rows repsresent grouping factors of interest. Note that the distribution might not be symmetrical, 
#'      in which case the \emph{se} is less informative.}
#' \item{CI_emp}{\code{data.frame} containing the (empirical) confidence intervals for the repeatabilities 
#'      estiamted based parametric bootstrapping. Each row represents a grouping factor of interest.}
#' \item{P}{\code{data.frame} with approximate p-values based on likelihood-ratio tests
#'      (first column) and permutation tests (second column). Each row represents a grouping factor 
#'      of interest.}
#' \item{R_boot}{Vector(s) of parametric bootstrap samples for \emph{R}. Each \code{list}
#'       element respesents a grouping factor.}
#' \item{R_permut}{Vector(s) of permutation samples for \emph{R}. Each \code{list}
#'       element represents a grouping factor.}
#' \item{LRT}{List of likelihoods for the full model and the reduced model(s), likelihood ratios \emph{D}, 
#'      \emph{P} value(s) and degrees of freedom for the likelihood-ratio test.} 
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
#' non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
#' 
#' @author Holger Schielzeth  (holger.schielzeth@@ebc.uu.se),
#'         Shinichi Nakagawa (shinichi.nakagawa@@otago.ac.nz) &
#'         Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'      
#' @seealso \link{rpt}
#' 
#' @examples  
#' 
#' data(BeetlesBody)
#' 
#' # Note: nboot and npermut are set to 5 for speed reasons. Use larger numbers
#' # for the real analysis.
#' 
#' # one random effect
#' rptGaussian(BodyL ~ (1|Population), grname="Population", 
#'                    data=BeetlesBody, nboot=5, npermut=5)
#' 
#' # two random effects
#' rptGaussian(BodyL ~ (1|Container) + (1|Population), grname=c("Container", "Population"), 
#'                    data=BeetlesBody, nboot=5, npermut=5)
#'                
#' 
#' 
#' @export
#' 

rptGaussian <- function(formula, grname, data, CI = 0.95, nboot = 1000, 
        npermut = 0, parallel = FALSE, ncores = NULL) {
        
        # missing values
        no_NA_vals <- stats::complete.cases(data[all.vars(formula)])
        if (sum(!no_NA_vals ) > 0 ){
                warning(paste0(sum(!no_NA_vals), " rows containing missing values were removed"))
                data <- data[no_NA_vals, ]
        } 
        
        mod <- lme4::lmer(formula, data = data)
        VarComps <- as.data.frame(lme4::VarCorr(mod))
        
        if (nboot == 1) {
                warning("nboot has to be greater than 1 to calculate a CI and has been set to 0")
                nboot <- 0
        }
        
        if (nboot < 0) nboot <- 0
        if (npermut < 1) npermut <- 1
        # point estimates of R
        R_pe <- function(formula, data, grname) {
                mod <- lme4::lmer(formula, data)
                VarComps <- lme4::VarCorr(mod)
                var_a <- as.numeric(VarComps[grname])
                names(var_a) <- grname
                var_p <- sum(as.numeric(VarComps)) + attr(VarComps, "sc")^2
                R <- var_a/var_p
                R <- as.data.frame(t(R))
                names(R) <- grname
                return(R)
        }
        
        R <- R_pe(formula, data, grname) # no bootstrap skipping at the moment
        
        # confidence interval estimation by parametric bootstrapping
        if (nboot > 0)  Ysim <- as.matrix(stats::simulate(mod, nsim = nboot))
        
        bootstr <- function(y, mod, formula, data, grname) {
                data[, names(stats::model.frame(mod))[1]] <- as.vector(y)
                R_pe(formula, data, grname)
        }
        if (nboot > 0 & parallel == TRUE) {
                if (is.null(ncores)) {
                        ncores <- parallel::detectCores() - 1
                        warning("No core number specified: detectCores() is used to detect the number of \n cores on the local machine")
                }
                # start cluster
                cl <- parallel::makeCluster(ncores)
                parallel::clusterExport(cl, "R_pe", envir=environment())
                R_boot <- unname(parallel::parApply(cl, Ysim, 2, bootstr, mod = mod, formula = formula, 
                        data = data, grname = grname))
                parallel::stopCluster(cl)
        }
        if (nboot > 0 & parallel == FALSE) {
                R_boot <- unname(apply(Ysim, 2, bootstr, mod = mod, formula = formula, data = data, 
                        grname = grname))
        }
        if (nboot == 0) {
                # R_boot <- matrix(rep(NA, length(grname)), nrow = length(grname))
                #                 R_boot <- list(structure(as.data.frame(matrix(rep(NA, 2*length(grname)), nrow = 2)),
                #                           names = grname, row.names = c("R_org", "R_link")))
                R_boot <- NA
        }
        
        # transform bootstrapping repeatabilities into vectors
        boot <- as.list(rep(NA, length(grname)))
        names(boot) <- grname
        # CI function
        calc_CI <- function(x) {
                out <- stats::quantile(x, c((1 - CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
        }
        
        if (length(R_boot) == 1) {
                # creating tables when R_boot = NA
                if (is.na(R_boot)) {
                       se <- NA
                       CI_emp <- calc_CI(NA)
                }
        } else  {
                for (i in 1:length(grname)) {
                        boot[[i]] <- unlist(lapply(R_boot, function(x) x[, grname[i]]))
                       
                }
                # CI 
                CI_emp <- as.data.frame(t(as.data.frame(lapply(boot, calc_CI))))
                # se
                se <- as.data.frame(t(as.data.frame(lapply(boot, stats::sd))))
                names(se) <- "se"
              
        }
        
        # significance test by permutation of residuals
        P_permut <- rep(NA, length(grname))
        
        # significance test by likelihood-ratio-test
        terms <- attr(terms(formula), "term.labels")
        randterms <- terms[which(regexpr(" | ", terms, perl = TRUE) > 0)]
        
        # no permutation test
        if (npermut == 1) {
                R_permut <- NA
                P_permut <- NA
        }
        
        # significance test by permutation of residuals
        # nperm argument just used for parallisation
        permut <- function(nperm, formula, data, mod_red, dep_var, grname, i) {
                y_perm <- stats::fitted(mod_red) + sample(stats::resid(mod_red))
                data_perm <- data
                data_perm[dep_var] <- y_perm
                out <- R_pe(formula, data_perm, grname[i])
                out
        }
        
        dep_var <- as.character(formula)[2]
        # one random effect, uses stats::lm()
        # multiple random effects, uses lmer()
 
        R_permut <- data.frame(matrix(rep(NA, length(grname) * npermut), nrow = length(grname)))
        P_permut <- rep(NA, length(grname))
        
        if (npermut > 1){
                for (i in 1:length(grname)) {
                        if (length(randterms) == 1) {
                                formula_red <- stats::update(formula, eval(paste(". ~ . ", paste("- (", randterms, ")"))))
                                mod_red <- stats::lm(formula_red, data = data)
                        } else if (length(randterms) > 1) {
                                formula_red <- stats::update(formula, eval(paste(". ~ . ", paste("- (1 | ", grname[i], 
                                        ")"))))
                                mod_red <- lme4::lmer(formula_red, data = data)
                        }
                        
                        if(parallel == TRUE) {
                                if (is.null(ncores)) {
                                        ncores <- parallel::detectCores()
                                        warning("No core number specified: detectCores() is used to detect the number of \n cores on the local machine")
                                }
                                # start cluster
                                cl <- parallel::makeCluster(ncores)
                                parallel::clusterExport(cl, "R_pe", envir=environment())
                                R_permut[i, ] <- c(R[i], as.numeric(unlist(parallel::parSapply(cl, 1:(npermut-1), permut, formula, data, mod_red, dep_var, grname, i))))
                                parallel::stopCluster(cl)
                                P_permut[i] <- sum(R_permut[i, ] >= unlist(R[i]))/npermut
                        } else if (parallel == FALSE) {
                                R_permut[i, ] <- c(R[i], as.numeric(unlist(replicate(npermut - 1, permut(formula=formula, data = data, 
                                        mod_red=mod_red, dep_var=dep_var, grname=grname, i=i), simplify = TRUE))))
                                P_permut[i] <- sum(R_permut[i, ] >= unlist(R[i]))/npermut
                        }
                }
        }
                
        row.names(R_permut) <- grname
        names(P_permut) <- grname
  
        
        ## likelihood-ratio-test
        LRT_mod <- as.numeric(stats::logLik(mod))
        LRT_df <- 1
        
        for (i in c("LRT_P", "LRT_D", "LRT_red")) assign(i, rep(NA, length(grname)))
        
        for (i in 1:length(grname)) {
                if (length(randterms) == 1) {
                        formula_red <- stats::update(formula, eval(paste(". ~ . ", paste("- (", randterms, ")"))))
                        LRT_red[i] <- as.numeric(stats::logLik(stats::lm(formula_red, data = data)))
                } else if (length(randterms) >= 1){
                        formula_red <- stats::update(formula, eval(paste(". ~ . ", paste("- (1 | ", grname[i], 
                        ")"))))
                        LRT_red[i] <- as.numeric(stats::logLik(lme4::lmer(formula = formula_red, data = data)))
                }
                LRT_D[i] <- as.numeric(-2 * (LRT_red[i] - LRT_mod))
                LRT_P[i] <- ifelse(LRT_D[i] <= 0, 1, stats::pchisq(LRT_D[i], 1, lower.tail = FALSE)/2)
                # LR <- as.numeric(-2*(logLik(lme4::lmer(stats::update(formula, eval(paste('. ~ . ',
                # paste('- (1 | ', grname[i], ')') ))), data=data))-logLik(mod))) P.LRT[i] <-
                # ifelse(LR<=0, 1, stats::pchisq(LR,1,lower.tail=FALSE)/2)
        }
        
        P <- cbind(LRT_P, P_permut)
    
        row.names(P) <- grname
        
        res <- list(call = match.call(), 
                datatype = "Gaussian", 
                CI = CI, 
                R = R, 
                se = se,
                CI_emp = CI_emp, 
                P = as.data.frame(P),
                R_boot = boot, 
                R_permut = lapply(as.data.frame(t(R_permut)), function(x) return(x)),
                LRT = list(LRT_mod = LRT_mod, LRT_red = LRT_red, LRT_D = LRT_D, LRT_df = LRT_df, 
                        LRT_P = LRT_P), 
                ngroups = unlist(lapply(data[grname], function(x) length(unique(x)))), 
                nobs = nrow(data), mod = mod)
        class(res) <- "rpt"
        return(res)
} 
