#' GLMM-based Repeatability Using REML
#' 
#' Calculates repeatability from a general linear mixed-effects models fitted by REML (restricted maximum likelihood).
#' @param formula Formula as used e.g. by \link{glmer}. The grouping factor(s) of
#'        interest needs to be included as a random effect, e.g. '(1|groups)'.
#'        Covariates and additional random effects can be included to estimate adjusted repeatabilities.
#' @param grname A character string or vector of character strings giving the
#'        name(s) of the grouping factor(s), for which the repeatability should
#'        be estimated. Spelling needs to match the random effect names as given in \code{formula}.
#' @param data A dataframe that contains the variables included in the \code{formula}
#'        and \code{grname} arguments.
#' @param link Link function. \code{log} and \code{sqrt} are allowed, defaults to \code{log}.
#' @param CI Width of the confidence interval (defaults to 0.95).
#' @param nboot Number of parametric bootstraps for interval estimation.
#'        Defaults to 1000. Larger numbers of permutations give a better
#'        asymtotic CI, but may be very time-consuming.
#' @param npermut Number of permutations used when calculating 
#'        asymptotic \emph{P} values (defaults to 1000). 
#' @param parallel If TRUE, bootstraps will be distributed. 
#' @param ncores Specify number of cores to use for parallelization. On default,
#'        all cores but one are used.
#' 
#' @return 
#' Returns an object of class rpt that is a a list with the following elements: 
#' \item{call}{function call}
#' \item{datatype}{Response distribution (here: 'Poisson').
#' \item{CI}{Width of the confidence interval.}
#' \item{R}{\code{data.frame} with point estimates for repeatabilities. Columns
#'      are groups of interest. Rows are original and link scale, in this order.}
#' \item{se}{\code{data.frame} with approximate standard errors (\emph{se}) for repeatabilities. Columns
#'      are groups of interest. Rows are original and link scale, in this order.
#'      Note that the distribution might not be symmetrical, in which case the \emph{se} is less informative.}
#' \item{CI_emp}{\code{list} of two elements containing the confidence intervals for repeatabilities 
#'      on the link and original scale, respectively. Within each list element, lower and upper CI
#'      are columns and each row is a grouping factor of interest.}
#' \item{P}{Approximate \emph{P} \code{data.frame} with p-values from a significance test based on likelihood-ratio
#'      in the first column and significance test based on permutation of residuals for 
#'      both the original and link scale in the second and third column. Each row is a grouping
#'      factors.}
#' \item{R_boot_link}{Parametric bootstrap samples for \emph{R} on the link scale. Each \code{list}
#'       element is a grouping factor.}
#' \item{R_boot_org}{Parametric bootstrap samples for \emph{R} on the original scale. Each \code{list}
#'       element is a grouping factor.}
#' \item{R_permut_link}{Permutation samples for \emph{R} on the link scale. Each \code{list}
#'       element is a grouping factor.}
#' \item{R_permut_org}{Permutation samples for \emph{R} on the original scale. Each \code{list}
#'       element is a grouping factor.}
#' \item{LRT}{List of Likelihood-ratios for the model and the reduced model(s), 
#'       and \emph{P} value(s) and degrees of freedom for the Likelihood-ratio test} 
#' \item{ngroups}{Number of groups.}
#' \item{nobs}{Number of observations.}
#' \item{overdisp}{Overdispersion parameter. Equals the variance in the observational factor random effect}
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
#' # repeatability estimation for tarsus length - a very high R
#' data(BodySize)
#' (rpt.BS <- rptGaussianNew(formula = Tarsus ~ 1 + (1|Sex) + (1|BirdID), grname = c('Sex', 'BirdID'), 
#'  data=BodySize, nboot=10, npermut=10))
#'  
#'  (rpt.BS <- rptGaussianNew(formula = Tarsus ~ 1 + (1|BirdID), grname = c('BirdID'), 
#'  data=BodySize, nboot=10, npermut=10))
#' @export
#' 

rptGaussianNew <- function(formula, grname, data, CI = 0.95, nboot = 1000, 
        npermut = 1000, parallel = FALSE, ncores = NULL) {
        
        # to do: missing values
        # no bootstrapping case
        
        mod <- lme4::lmer(formula, data = data)
        VarComps <- as.data.frame(lme4::VarCorr(mod))

        # check if all variance components are 0
#         if (sum(VarComps$vcov[-"Residual"]!=0) == 0) {
#                 nboot <- 0
#                 npermut <- 0
#                 warning("all variance components are 0, bootstrapping and permutation skipped")
#         }
        
        
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
                data[, names(model.frame(mod))[1]] <- as.vector(y)
                R_pe(formula, data, grname)
        }
        if (nboot > 0 & parallel == TRUE) {
                if (is.null(ncores)) {
                        ncores <- parallel::detectCores() - 1
                        warning("No core number specified: detectCores() is used to detect the number of \n cores on the local machine")
                }
                # start cluster
                cl <- parallel::makeCluster(ncores)
                parallel::clusterExport(cl, "R_pe")
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
                out <- quantile(x, c((1 - CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
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
                se <- as.data.frame(t(as.data.frame(lapply(boot, sd))))
                names(se) <- "se"
              
        }
        
        # significance test by permutation of residuals
        P_permut <- rep(NA, length(grname))
        
        # significance test by likelihood-ratio-test
        terms <- attr(terms(formula), "term.labels")
        randterms <- terms[which(regexpr(" | ", terms, perl = TRUE) > 0)]
        
        # no permutation test
        if (npermut == 1) {
                R_permut <- R
                P_permut <- NA
        }
        
        # significance test by permutation of residuals
        # nperm argument just used for parallisation
        permut <- function(nperm, formula, data, mod_red, dep_var, grname, i) {
                y_perm <- fitted(mod_red) + sample(resid(mod_red))
                data_perm <- data
                data_perm[dep_var] <- y_perm
                out <- R_pe(formula, data_perm, grname[i])
                out
        }
        
        # multiple random effects, uses lmer()
        dep_var <- as.character(formula)[2]
        # one random effect, uses lm()
        # multiple random effects, uses lmer()
 
        R_permut <- data.frame(matrix(rep(NA, length(grname) * npermut), nrow = length(grname)))
        P_permut <- rep(NA, length(grname))
        for (i in 1:length(grname)) {
                if (length(randterms) == 1) {
                        formula_red <- update(formula, eval(paste(". ~ . ", paste("- (", randterms, ")"))))
                        mod_red <- lm(formula_red, data = data)
                } else if (length(randterms) > 1) {
                        formula_red <- update(formula, eval(paste(". ~ . ", paste("- (1 | ", grname[i], 
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
                        parallel::clusterExport(cl, "R_pe")
                        R_permut[i, ] <- c(R[i], as.numeric(unlist(parallel::parSapply(cl, 1:(npermut-1), permut, formula, data, mod_red, dep_var, grname, i))))
                        parallel::stopCluster(cl)
                        P_permut[i] <- sum(R_permut[i, ] >= unlist(R[i]))/npermut
                } else if (parallel == FALSE) {
                        R_permut[i, ] <- c(R[i], as.numeric(unlist(replicate(npermut - 1, permut(formula=formula, data = data, 
                                mod_red=mod_red, dep_var=dep_var, grname=grname, i=i), simplify = TRUE))))
                        P_permut[i] <- sum(R_permut[i, ] >= unlist(R[i]))/npermut
                }
        }
                
        row.names(R_permut) <- grname
        names(P_permut) <- grname
  
        
        ## likelihood-ratio-test
        LRT_mod <- as.numeric(logLik(mod))
        LRT_df <- 1
        
        for (i in c("LRT_P", "LRT_D", "LRT_red")) assign(i, rep(NA, length(grname)))
        
        for (i in 1:length(grname)) {
                if (length(randterms) == 1) {
                        formula_red <- update(formula, eval(paste(". ~ . ", paste("- (", randterms, ")"))))
                        LRT_red[i] <- as.numeric(logLik(lm(formula_red, data = data)))
                } else if (length(randterms) >= 1){
                        formula_red <- update(formula, eval(paste(". ~ . ", paste("- (1 | ", grname[i], 
                        ")"))))
                        LRT_red[i] <- as.numeric(logLik(lme4::lmer(formula = formula_red, data = data)))
                }
                LRT_D[i] <- as.numeric(-2 * (LRT_red[i] - LRT_mod))
                LRT_P[i] <- ifelse(LRT_D[i] <= 0, 1, pchisq(LRT_D[i], 1, lower.tail = FALSE)/2)
                # LR <- as.numeric(-2*(logLik(lme4::lmer(update(formula, eval(paste('. ~ . ',
                # paste('- (1 | ', grname[i], ')') ))), data=data))-logLik(mod))) P.LRT[i] <-
                # ifelse(LR<=0, 1, pchisq(LR,1,lower.tail=FALSE)/2)
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
