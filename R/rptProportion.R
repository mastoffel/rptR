#' GLMM-based Repeatability Using REML for Proportion data
#' 
#' Calculates repeatability from a general linear mixed-effects models fitted by REML (restricted maximum likelihood).
#' 
#' Calculates repeatability from a generalized linear mixed-effects models fitted by restricted maximum likelihood (REML).
#' @param formula Formula as used e.g. by \link{lmer}. In the case of Porportion data, the left-hand
#'        side of the formula represented to vectors with counts of successes and failures, 
#'        respectively, connected by \code{cbind}. The grouping factor(s) of
#'        interest needs to be included as a random effect, e.g. '(1|groups)'.
#'        Covariates and additional random effects can be included to estimate adjusted 
#'        repeatabilities.
#' @param grname A character string or vector of character strings giving the
#'        name(s) of the grouping factor(s), for which the repeatability should
#'        be estimated. Spelling needs to match the random effect names as given in \code{formula} 
#'        and terms have to be set in quotation marks.
#' @param data A dataframe that contains the variables included in the \code{formula}
#'        and \code{grname} arguments.
#' @param link Link function. \code{logit} and \code{probit} are allowed, defaults to \code{logit}.
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
#' \item{datatype}{Response distribution (here: 'Proportion')}.
#' \item{CI}{Width of the confidence interval.}
#' \item{R}{\code{data.frame} with point estimates for repeatabilities. Columns
#'      represent grouping factors of interest. Rows show original and link scale repeatabilites 
#'      (in this order).}
#' \item{se}{\code{data.frame} with approximate standard errors (\emph{se}) for repeatabilities. Columns
#'      are groups of interest. Rows are original and link scale, in this order.
#'      Note that the distribution might not be symmetrical, in which case the \emph{se} is less informative.}
#' \item{CI_emp}{\code{list} of two elements containing the confidence intervals for repeatabilities 
#'      on the link and original scale, respectively. Within each list element, lower and upper CI
#'      are columns and each row for each grouping factor of interest.}
#' \item{P}{Approximate \emph{P} \code{data.frame} with p-values from a significance test based on likelihood-ratio
#'      in the first column and significance test based on permutation of residuals for 
#'      both the original and link scale in the second and third column. Each row represents a grouping
#'      factor of interest.}
#' \item{R_boot_link}{Parametric bootstrap samples for \emph{R} on the link scale. Each \code{list}
#'       element is a grouping factor.}
#' \item{R_boot_org}{Parametric bootstrap samples for \emph{R} on the original scale. Each \code{list}
#'       element is a grouping factor.}
#' \item{R_permut_link}{Permutation samples for \emph{R} on the link scale. Each \code{list}
#'       element is a grouping factor.}
#' \item{R_permut_org}{Permutation samples for \emph{R} on the original scale. Each \code{list}
#'       element is a grouping factor.}
#' \item{LRT}{List of likelihoods for the full model and the reduced model(s), likelihood ratios \emph{D}, 
#'      \emph{P} value(s) and degrees of freedom for the likelihood-ratio test.} 
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
#'  
#' @export
#' 

rptProportion <- function(formula, grname, data, link = c("logit", "probit"), CI = 0.95, nboot = 1000, 
        npermut = 1000, parallel = FALSE, ncores = NULL) {
        
        # missing values
        no_NA_vals <- stats::complete.cases(data[all.vars(formula)])
        if (sum(!no_NA_vals ) > 0 ){
                warning(paste0(sum(!no_NA_vals), " rows containing missing values were removed"))
                data <- data[no_NA_vals, ]
        } 
        
        # link
        if (length(link) > 1) link <- "logit"
        if (!(link %in% c("logit", "probit"))) stop("Link function has to be 'logit' or 'probit'")
        # observational level random effect
        obsid <- factor(1:nrow(data))
        data <- cbind(data, obsid)
        
        formula <- stats::update(formula,  ~ . + (1|obsid))
        mod <- lme4::glmer(formula, data = data, family = stats::binomial(link = link))
        VarComps <- as.data.frame(lme4::VarCorr(mod))
        obsind_id <- which(VarComps[["grp"]] == "obsid")
        overdisp <- VarComps$vcov[obsind_id]
        
        # check if all variance components are 0
        if (sum(VarComps$vcov[-obsind_id]!=0) == 0) {
                nboot <- 0
                npermut <- 0
                warning("all variance components are 0, bootstrapping and permutation skipped")
        }
        
        if (nboot == 1) {
                warning("nboot has to be greater than 1 to calculate a CI and has been set to 0")
                nboot <- 0
        }
        
        if (nboot < 0) nboot <- 0
        if (npermut < 1) npermut <- 1
        e1 <- environment()
        # point estimates of R
        R_pe <- function(formula, data, grname, peYN = FALSE) {
                
                mod <- suppressWarnings(lme4::glmer(formula = formula, data = data, family = stats::binomial(link = link)))
                # random effect variance data.frame
                VarComps <- as.data.frame(lme4::VarCorr(mod))
                # groups random effect variances
                var_a <- VarComps[VarComps$grp %in% grname, "vcov"]
                names(var_a) <- grname
                # olre variance
                var_e <- VarComps[VarComps$grp %in% "obsid", "vcov"]
                # intercept on link scale
                beta0 <- unname(lme4::fixef(mod)[1])
  
                if (link == "logit") {
                        R_link <- var_a/(var_a + var_e + (pi^2)/3)
                        P <- exp(beta0) / (1 + exp(beta0))
                        R_org <- ( (var_a * P^2) / ((1 + exp(beta0))^2)) / (((var_a + var_e) * P^2) / ((1 + exp(beta0))^2) + (P * (1-P)))
                }
                if (link == "probit") {
                        R_link <- var_a/(var_a + var_e + 1)
                        R_org <- NA
                        
                }
                # check whether that works for any number of var
                R <- as.data.frame(rbind(R_org, R_link))
                return(R)
                }
        
        R <- R_pe(formula, data, grname, peYN = FALSE)
        
        # confidence interval estimation by parametric bootstrapping
        if (nboot > 0){
                Ysim <- stats::simulate(mod, nsim = nboot) # every column contains a list with successes and failures
        }  
        
        bootstr <- function(y, mod, formula, data, grname) {
                data[, colnames(y)] <- y
                R_pe(formula, data, grname)
        }
        
        # to do: preallocate R_boot
        if (nboot > 0 & parallel == TRUE) {
                if (is.null(ncores)) {
                        ncores <- parallel::detectCores() - 1
                        warning("No core number specified: detectCores() is used to detect the number of \n cores on the local machine")
                }
                # start cluster
                cl <- parallel::makeCluster(ncores)
                parallel::clusterExport(cl, "R_pe")
                R_boot <- unname(parallel::parLapply(cl, Ysim, bootstr, mod = mod, formula = formula, 
                        data = data, grname = grname))
                parallel::stopCluster(cl)
        }
        if (nboot > 0 & parallel == FALSE) {
                R_boot <- unname(lapply(Ysim, bootstr, mod = mod, formula = formula, data = data, 
                        grname = grname))
        }
        if (nboot == 0) {
                # R_boot <- matrix(rep(NA, length(grname)), nrow = length(grname))
                R_boot <- NA
        }
        
        # transform bootstrapping repeatabilities into vectors
        boot_org <- as.list(rep(NA, length(grname)))
        boot_link <- as.list(rep(NA, length(grname)))
        if (length(R_boot) == 1) {
                # creating tables when R_boot = NA
                if (is.na(R_boot)) {
                        # for(i in c("CI_org", "CI_link", "se_org", "se_link")) assign(i, NA, envir = e1)
                        for(i in c("se_org", "se_link")){
                                assign(i, structure(data.frame(matrix(NA, 
                                        nrow = length(grname))), row.names = grname, names = i), 
                                        envir = e1)   
                        }
                        for(i in c("CI_org", "CI_link")){
                                assign(i, structure(data.frame(matrix(NA, 
                                        nrow = length(grname), ncol = 2)), row.names = grname), 
                                        envir = e1)   
                        }
                        
                }
        }  else {
                for (i in 1:length(grname)) {
                        boot_org[[i]] <- unlist(lapply(R_boot, function(x) x["R_org", grname[i]]))
                        boot_link[[i]] <- unlist(lapply(R_boot, function(x) x["R_link", grname[i]]))
                }
                names(boot_org) <- grname
                names(boot_link) <- grname
                
                calc_CI <- function(x) {
                        out <- stats::quantile(x, c((1 - CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
                }
                
                # CI into data.frame and transpose to have grname in rows
                CI_org <- as.data.frame(t(as.data.frame(lapply(boot_org, calc_CI))))
                CI_link <- as.data.frame(t(as.data.frame(lapply(boot_link, calc_CI))))
                
                # se
                se_org <- as.data.frame(t(as.data.frame(lapply(boot_org, stats::sd))))
                se_link <- as.data.frame(t(as.data.frame(lapply(boot_link, stats::sd))))
                names(se_org) <- "se_org"
                names(se_link) <- "se_link"
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
        
        # response matrix for permutation test
        dep_var_expr <- as.character(formula)[2]
        dep_var <- as.data.frame(with(data, eval(parse(text=dep_var_expr))))
        
        # significance test by permutation of residuals
        # nperm argument just used for parallisation
        
        if (link == "logit") {
                trans_fun <- stats::qlogis     # VGAM::logit
                inv_fun <- stats::plogis
        }
        if (link == "probit") {
                trans_fun <- stats::qnorm            # VGAM::probit
                inv_fun <- stats::pnorm
        }
        
        permut <- function(nperm, formula, mod, dep_var, grname, data) {
                # for binom it will be logit 
                 y_perm <- stats::rbinom(nrow(data), rowSums(dep_var), 
                                  prob = inv_fun((trans_fun(stats::fitted(mod)) + 
                                  sample(stats::resid(mod)))))
                # y_perm <- rbinom(nrow(data), 1, prob = (predict(mod, type = "response") + sample(resid(mod))))
                data_perm <- data
                data_perm[names(dep_var)[1]] <- y_perm
                data_perm[names(dep_var)[2]] <- rowSums(dep_var) - y_perm
                out <- R_pe(formula, data_perm, grname)
                out
        }
        
        # R_permut <- matrix(rep(NA, length(grname) * npermut), nrow = length(grname))
        P_permut <- structure(data.frame(matrix(NA, nrow = 2, ncol = length(grname)),
                row.names = c("P_permut_org", "P_permut_link")), names = grname)
        
        # for likelihood ratio and permutation test
        terms <- attr(terms(formula), "term.labels")
        randterms <- terms[which(regexpr(" | ", terms, perl = TRUE) > 0)]
        
        if (npermut > 1){
                for (i in 1:length(grname)) {
                        if (length(randterms) > 1) {
                                formula_red <- stats::update(formula, eval(paste(". ~ . ", paste("- (1 | ", grname[i], 
                                        ")"))))
                                mod_red <- lme4::glmer(formula_red, data = data, family = stats::binomial(link = link))
                        } else if (length(randterms) == 1) {
                                formula_red <- stats::update(formula, eval(paste(". ~ . ", paste("- (", randterms, ")"))))
                                mod_red <- stats::glm(formula_red, data = data, family = stats::binomial(link = link))
                        }
                        
                if(parallel == TRUE) {
                        if (is.null(ncores)) {
                                ncores <- parallel::detectCores()
                                warning("No core number specified: detectCores() is used to detect the number of \n cores on the local machine")
                        }
                # start cluster
                cl <- parallel::makeCluster(ncores)
                parallel::clusterExport(cl, "R_pe")
                R_permut <- parallel::parLapply(cl, 1:(npermut-1), permut, formula=formula, 
                        mod=mod_red, dep_var=dep_var, grname=grname, data = data)
                parallel::stopCluster(cl)
                
                } else if (parallel == FALSE) {
                        R_permut <- lapply(1:(npermut - 1), permut, formula, mod_red, dep_var, grname, data)
                }
                # adding empirical rpt 
                R_permut <- c(list(R), R_permut)
                }
        }
        
        # equal to boot
        permut_org <- as.list(rep(NA, length(grname)))
        permut_link <- as.list(rep(NA, length(grname)))
        
        if (!(length(R_permut) == 1)){
                        for (i in 1:length(grname)) {
                                permut_org[[i]] <- unlist(lapply(R_permut, function(x) x["R_org", grname[i]]))
                                permut_link[[i]] <- unlist(lapply(R_permut, function(x) x["R_link", grname[i]]))
                        }
                        names(permut_org) <- grname
                        names(permut_link) <- grname
        }
        
        
        P_permut["P_permut_org", ] <- unlist(lapply(permut_org, function(x) sum(x >= x[1])))/npermut
        P_permut["P_permut_link", ] <- unlist(lapply(permut_link, function(x) sum(x >= x[1])))/npermut
        names(P_permut) <- names(permut_link)
        
        
        ## likelihood-ratio-test
        LRT_mod <- as.numeric(stats::logLik(mod))
        LRT_df <- 1
        #         if (length(randterms) == 1) {
        #                 formula_red <- update(formula, eval(paste(". ~ . ", paste("- (", randterms, ")"))))
        #                 LRT.red <- as.numeric(logLik(lm(formula_red, data = data)))
        #                 LRT.D <- as.numeric(-2 * (LRT.red - LRT.mod))
        #                 LRT.P <- ifelse(LRT.D <= 0, LRT.df, pchisq(LRT.D, 1, lower.tail = FALSE)/2)
        #                 # LR <- as.numeric(-2*(logLik(lm(update(formula, eval(paste('. ~ . ', paste('- (',
        #                 # randterms, ')') ))), data=data))-logLik(mod))) P.LRT <- ifelse(LR<=0, 1,
        #                 # pchisq(LR,1,lower.tail=FALSE)/2)
        #         }
        
        
        for (i in c("LRT_P", "LRT_D", "LRT_red")) assign(i, rep(NA, length(grname)))
        
        for (i in 1:length(grname)) {
                formula_red <- stats::update(formula, eval(paste(". ~ . ", paste("- (1 | ", grname[i], 
                        ")"))))
                LRT_red[i] <- as.numeric(stats::logLik(lme4::glmer(formula = formula_red, data = data, 
                        family = stats::binomial(link = link))))
                LRT_D[i] <- as.numeric(-2 * (LRT_red[i] - LRT_mod))
                LRT_P[i] <- ifelse(LRT_D[i] <= 0, 1, stats::pchisq(LRT_D[i], 1, lower.tail = FALSE)/2)
                # LR <- as.numeric(-2*(logLik(lme4::lmer(update(formula, eval(paste('. ~ . ',
                # paste('- (1 | ', grname[i], ')') ))), data=data))-logLik(mod))) P.LRT[i] <-
                # ifelse(LR<=0, 1, pchisq(LR,1,lower.tail=FALSE)/2)
        }
        
        P <- cbind(LRT_P, t(P_permut))
        row.names(P) <- grname
        
        res <- list(call = match.call(), 
                datatype = "Binary", 
                link = link,
                CI = CI, 
                R = R, 
                se = t(cbind(se_org,se_link)), 
                CI_emp = list(CI_org = CI_org, CI_link = CI_link), 
                P = as.data.frame(P),
                R_boot_link = boot_link, 
                R_boot_org = boot_org,
                R_permut_link = permut_link, 
                R_permut_org = permut_org,
                LRT = list(LRT_mod = LRT_mod, LRT_red = LRT_red, LRT_D = LRT_D, LRT_df = LRT_df, 
                        LRT_P = LRT_P), 
                ngroups = unlist(lapply(data[grname], function(x) length(unique(x)))), 
                nobs = nrow(data), overdisp = overdisp, mod = mod)
        class(res) <- "rpt"
        return(res)
} 