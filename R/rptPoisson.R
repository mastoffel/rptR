#' GLMM-based Repeatability Estimation for Poisson-distributed Data
#' 
#' Estimates repeatability from a generalized linear mixed-effects models fitted by restricted maximum likelihood (REML).
#' @inheritParams rpt
#' @param link Link function. \code{logit} and \code{probit} are allowed, defaults to \code{logit}.
#' 
#' 
#' @return 
#' Returns an object of class \code{rpt} that is a a list with the following elements: 
#' \item{call}{Function call}
#' \item{datatype}{Response distribution (here: 'Poisson').}
#' \item{CI}{Coverage of the confidence interval as specified by the \code{CI} argument.}
#' \item{R}{\code{data.frame} with point estimates for repeatabilities. Columns
#'      represent grouping factors of interest. Rows show original and link scale repeatabilites 
#'      (in this order).}
#' \item{se}{\code{data.frame} with approximate standard errors (\emph{se}) for repeatabilities. Columns
#'      are groups of interest. Rows are original and link scale (in this order).
#'      Note that the distribution might not be symmetrical, in which case the \emph{se} is less informative.}
#' \item{CI_emp}{\code{list} of two elements containing the confidence intervals for repeatabilities 
#'      on the link and original scale, respectively. Within each list element, lower and upper CI
#'      are columns and each row for each grouping factor of interest.}
#' \item{P}{\code{data.frame} with p-values from a significance test based on likelihood-ratios
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
#' \item{LRT}{List of two elements. \emph{LRT_mod} is the likelihood for the full model and (2) \emph{LRT_table} is a data.frame 
#'      for the reduced model(s) including columns for the likelihood \emph{logl_red}, the likelihood ratio(s) \emph{LR_D}, 
#'      p-value(s)\emph{LR_P} and degrees of freedom for the likelihood-ratio test(s) \emph{LR_df}.} 
#' \item{ngroups}{Number of groups for each grouping level.}
#' \item{nobs}{Number of observations.}
#' \item{mod}{Fitted model.}
#' \item{ratio}{Boolean. TRUE, if ratios have been estimated, FALSE, if variances have been estimated}
#' \item{adjusted}{Boolean. TRUE, if estimates are adjusted}
#' \item{all_warnings}{\code{list} with two elements. 'warnings_boot' and 'warnings_permut' contain
#'      warnings from the lme4 model fitting of bootstrap and permutation samples, respectively.}
#'
#' @details 
#' 
#' see details section of \code{\link{rpt}} for details on parametric bootstrapping,
#' permutation and likelihood-ratio tests.
#' 
#' @references 
#' Carrasco, J. L. & Jover, L.  (2003) \emph{Estimating the generalized 
#' concordance correlation coefficient through variance components}. Biometrics 59: 849-858.
#'
#' Faraway, J. J. (2006) \emph{Extending the linear model with R}. Boca Raton, FL, Chapman & Hall/CRC.
#' 
#' Nakagawa, S. & Schielzeth, H. (2010) \emph{Repeatability for Gaussian and 
#' non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
#' 
#' @author Holger Schielzeth  (holger.schielzeth@@uni-jena.de),
#'         Shinichi Nakagawa (s.nakagawa@unsw.edu.au) &
#'         Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'      
#' @seealso \link{rpt}
#' 
#' 
#' @examples 
#' # load data
#' data(BeetlesFemale)
#' 
#' # Note: nboot and npermut are set to 0 for speed reasons. 
#' 
#' # estimating adjusted repeatabilities for two random effects
#' rptPoisson(Egg ~ Treatment + (1|Container) + (1|Population), 
#'                    grname=c("Container", "Population"), 
#'                    data = BeetlesFemale, nboot=0, npermut=0)
#' 
#' # unadjusted repeatabilities with  fixed effects and 
#' # estimation of the fixed effect variance
#' rptPoisson(Egg ~ Treatment + (1|Container) + (1|Population), 
#'                    grname=c("Container", "Population", "Fixed"), 
#'                    data=BeetlesFemale, nboot=0, npermut=0, adjusted=FALSE)
#'                 
#'                
#' # variance estimation of random effects, residual and overdispersion 
#' rptPoisson(formula = Egg ~ Treatment + (1|Container) + (1|Population) , 
#'                    grname=c("Container","Population","Residual", "Overdispersion"), 
#'                    data = BeetlesFemale, nboot=0, npermut=0, ratio = FALSE)
#'                    
#'                    
#' 
#' 
#' @export
#' 

rptPoisson <- function(formula, grname, data, link = c("log", "sqrt"), CI = 0.95, nboot = 1000, 
        npermut = 0, parallel = FALSE, ncores = NULL, ratio = TRUE, adjusted = TRUE, expect="meanobs",
        rptOutput = NULL, update = FALSE) {
        
        # missing values
        no_NA_vals <- stats::complete.cases(data[all.vars(formula)])
        if (sum(!no_NA_vals ) > 0 ){
                warning(paste0(sum(!no_NA_vals), " rows containing missing values were removed"))
                data <- data[no_NA_vals, ]
        } 
        
        # check whether grnames just contain "Residual" or "Overdispersion"
        if (!any((grname != "Residual") & (grname != "Fixed"))) stop("Specify at least one grouping factor in grname")
        
        # link
        if (length(link) > 1) link <- "log" 
        if (!(link %in% c("log", "sqrt"))) stop("Link function has to be 'log' or 'sqrt'")
        
        # check whether expect is either "meanobs" or "latent"
        if (link == "log" & (expect != "meanobs" & expect != "latent")) stop("The argument expect has to be either 'meanobs' (the default) or 'latent'")
        
        # observational level random effect
        Overdispersion <- factor(1:nrow(data))
        data <- cbind(data, Overdispersion)
        formula <- stats::update(formula,  ~ . + (1|Overdispersion))
        mod <- lme4::glmer(formula, data = data, family = stats::poisson(link = link))
        
        if (nboot == 1) {
                warning("nboot has to be greater than 1 to calculate a CI and has been set to 0")
                nboot <- 0
        }
        if (nboot < 0) nboot <- 0
        if (npermut < 1) npermut <- 1
        e1 <- environment()
        
        # save the original grname
        grname_org <- grname
        output_resid <- FALSE
        output_fixed <- FALSE
        
        # check whether Residual, Overdispersion or Fixed is selected and if so, remove it
        # from grname vector
        for (component in c("Residual", "Fixed")) {
                if (any(grname == component)){
                        grname <- grname[-which(grname == component)]
                        if (component == "Residual") output_resid <- TRUE
                        if (component == "Fixed") output_fixed <- TRUE
                }
        }
        
        # point estimates of R
        R_pe <- function(formula, data, grname, peYN = FALSE) {
                
                # mod <- suppressWarnings(lme4::glmer(formula = formula, data = data, family = stats::poisson(link = link)))
                mod <- lme4::glmer(formula = formula, data = data, family = stats::poisson(link = link))
                
                # random effect variance data.frame
                VarComps <- as.data.frame(lme4::VarCorr(mod))
                rownames(VarComps) = VarComps$grp
                
                # groups random effect variances
                var_a <- VarComps[grname, "vcov"]
                names(var_a) <- grname
                
                # intercept on link scale
                beta0 <- unname(lme4::fixef(mod)[1])
                
                # Fixed effect variance
                var_f <- stats::var(stats::predict(mod, re.form=NA))
                
                # Distribution-specific and residual variance
                if (link == "sqrt") {
                        estdv_link = 0.25
                        var_r <- VarComps["Overdispersion", "vcov"] + estdv_link
                }
                if (link == "log") {
                        if(expect=="meanobs") EY <- mean(mod@resp$y, na.rm=TRUE)
                        if(expect=="latent") EY <- exp(beta0 + (sum(VarComps[,"vcov"]) + var_f)/2)
                        estdv_link = log(1/EY+1)
                        var_r <- VarComps["Overdispersion", "vcov"] + estdv_link
                }
                
                if (ratio == FALSE) {
                        R_link <- var_a
                        R_org <- NA
                        R <- as.data.frame(rbind(R_org, R_link))
                        # if residual is selected, add residual variation to the output
                        if (output_resid){
                                R[,"Residual"] <- c(NA,var_r)  # add NA for R_org
                        } 
                        if (output_fixed){
                                R[,"Fixed"] <- c(NA, var_f)  # add NA for R_org
                        } 
                        return(R)
                }
                
                # Repeatability
                if (ratio == TRUE) {
                        if (link == "sqrt") {
                                # link scale
                                var_p_link <- sum(VarComps[,"vcov"]) + estdv_link
                                if(!adjusted) var_p_link <- var_p_link + var_f
                                R_link <- var_a / var_p_link
                                R_r <- var_r / var_p_link
                                R_f_link <- var_f / var_p_link
                                # origial scale
                                R_org <- NA
                                R_f_org <- NA
                        }
                        if (link == "log") {
                                # link scale
                                var_p_link <- sum(VarComps[,"vcov"]) +  estdv_link
                                if(!adjusted) var_p_link <- var_p_link + var_f
                                R_link = var_a / var_p_link
                                R_r <- var_r / var_p_link
                                R_f_link <- var_f / var_p_link
                                # origial scale
                                if( adjusted) var_p_org <- EY * (exp(sum(VarComps[,"vcov"])) - 1) + 1
                                if(!adjusted) var_p_org <- EY * (exp(sum(VarComps[,"vcov"] + var_f)) - 1) + 1
                                R_org <- EY * (exp(var_a) - 1)/ var_p_org
                                R_f_org <- EY * (exp(var_f) - 1)/ var_p_org
                        }
                        # check whether that works for any number of var
                        R <- as.data.frame(rbind(R_org, R_link))
                        
                        # check whether to give out non-repeatability and overdispersion repeatability
                        if (output_resid){
                                R[,"Residual"] <- c(NA,R_r) # add NA for R_org
                        }
                        if (output_fixed){
                                R[,"Fixed"] <- c(R_f_org, R_f_link) 
                        }
                        
                        return(R)
                }
        }
        
        
        R <- R_pe(formula, data, grname, peYN = FALSE) # no bootstrap skipping at the moment
        
        # confidence interval estimation by parametric bootstrapping
        
        # simulation of data.frame with responses
        if (nboot > 0)  Ysim <- as.data.frame(stats::simulate(mod, nsim = nboot))
        # main bootstrap function
        bootstr <- function(y, mod, formula, data, grname) {
                data[, names(stats::model.frame(mod))[1]] <- as.vector(y)
                R_pe(formula, data, grname)
        }
        
        # run all bootstraps
        bootstraps <- bootstrap_nongaussian(bootstr, R_pe, formula, data, Ysim, mod, grname, 
                grname_org, nboot, parallel, ncores, CI, rptOutput, update)
        
        # load everything (elegant solution)
        # list2env(bootstraps, envir = e1)
        
        # load everything (bad solution to assure global binding and satisfy cran check) 
        se_org <- bootstraps$se_org
        se_link <- bootstraps$se_link
        CI_org <- bootstraps$CI_org
        CI_link <- bootstraps$CI_link
        boot_link <- bootstraps$boot_link
        boot_org <- bootstraps$boot_org
        warnings_boot <- bootstraps$warnings_boot
        
        
        
        ### significance test by permutation of residuals ###
        
        # response variable
        dep_var <- as.character(formula)[2]
        
        #  main permutation function
        permut <- function(nperm, formula, mod_red, dep_var, grname, data) {
                if (link == "sqrt") {
                        y_perm <- stats::rpois(nrow(data), 
                                (stats::predict(mod_red, type="link") + sample(stats::resid(mod_red)))^2)
                }
                if (link == "log") {
                        y_perm <- stats::rpois(nrow(data), 
                                exp(stats::predict(mod_red, type="link") + sample(stats::resid(mod_red))))
                }
                data_perm <- data
                data_perm[dep_var] <- y_perm
                out <- R_pe(formula, data_perm, grname)
                out
        }
        
        family <- "poisson"
        permutations <- permut_nongaussian(permut, R_pe, formula, data, dep_var, 
                grname, npermut, parallel, ncores, link, family, R, rptOutput, update)
        
        P_permut <- permutations$P_permut
        permut_org <- permutations$permut_org
        permut_link <- permutations$permut_link
        warnings_permut <- permutations$warnings_permut
        
        
        
        ### likelihood-ratio-test ###
        LRTs <- LRT_nongaussian(formula, data, grname, mod, link, family)
        
        LRT_mod <- LRTs$LRT_mod
        LRT_table <- LRTs$LRT_table
        LRT_P <- LRT_table$LRT_P
        P <- cbind(LRT_P, t(P_permut))
        row.names(P) <- grname
        
        
        # add Residual = NA for S3 functions to work
        for (component in c("Residual", "Fixed")) {
                if(any(grname_org == component)){
                        # grname <- grname_org
                        P <- rbind(P, as.numeric(NA))
                        row.names(P)[nrow(P)] <- component
                        permut_link[component] <- NA
                        permut_org[component] <- NA
                        LRT_table <- rbind(LRT_table, as.numeric(NA))
                        row.names(LRT_table)[nrow(LRT_table)] <- component
                }
        }
        
        # delete overdispersion from ngroups
        ngroups <-  unlist(lapply(data[grname], function(x) length(unique(x))))
        ngroups <- ngroups[!names(ngroups) == "Overdispersion"]
        
        res <- list(call = match.call(), 
                datatype = "Poisson", 
                link = link,
                CI = CI, 
                R = R, 
                se = as.data.frame(t(cbind(se_org,se_link))),
                CI_emp = list(CI_org = CI_org, CI_link = CI_link), 
                P = as.data.frame(P),
                R_boot_link = boot_link, 
                R_boot_org = boot_org,
                R_permut_link = permut_link, 
                R_permut_org = permut_org,
                LRT = list(LRT_mod = LRT_mod, LRT_table = LRT_table), 
                ngroups =ngroups, 
                nobs = nrow(data), mod = mod, ratio = ratio, adjusted = adjusted,
                all_warnings = list(warnings_boot = warnings_boot, warnings_permut = warnings_permut))
        
        class(res) <- "rpt"
        return(res)
} 