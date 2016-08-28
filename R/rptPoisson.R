#' GLMM-based Repeatability Estimation for Possoin-distributed Data
#' 
#' Estimates repeatability from a generalized linear mixed-effects models fitted by restricted maximum likelihood (REML).
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
#' @param link Link function. \code{logit} and \code{probit} are allowed, defaults to \code{logit}.
#' @param CI Width of the required confidence interval (defaults to 0.95).
#' @param nboot Number of parametric bootstraps for interval estimation 
#'        (defaults to 1000). Larger numbers of bootstraps give a better
#'        asymtotic CI, but may be time-consuming. Bootstrapping can be switch off by setting 
#'        \code{nboot = 0}.
#' @param npermut Number of permutations used when calculating asymptotic p-values 
#'        (defaults to 0). Larger numbers of permutations give a better
#'        asymtotic p-values, but may be time-consuming (in particular when multiple grouping factors
#'        are specified). Permutaton tests can be switch off by setting \code{npermut = 0}. 
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
#' \item{LRT}{List of likelihoods for the full model and the reduced model(s), likelihood ratios \emph{D}, 
#'      p-value(s) and degrees of freedom for the likelihood-ratio test.} 
#' \item{ngroups}{Number of groups for each grouping level.}
#' \item{nobs}{Number of observations.}
#' \item{mod}{Fitted model.}
#' \item{all_warnings}{\code{list} with two elements. 'warnings_boot' and 'warnings_permut' contain
#'      warnings from the lme4 model fitting of bootstrap and permutation samples, respectively.}
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
#' # Note: nboot and npermut are set to 3 for speed reasons. Use larger numbers
#' # for the real analysis.
#' 
#' # estimating adjusted repeatabilities for two random effects
#' rptPoisson(Egg ~ Treatment + (1|Container) + (1|Population), 
#'                    grname=c("Container", "Population"), 
#'                    data = BeetlesFemale, nboot=3, npermut=3)
#' 
#' # unadjusted repeatabilities with  fixed effects and 
#' # estimation of the fixed effect variance
#' rptPoisson(Egg ~ Treatment + (1|Container) + (1|Population), 
#'                    grname=c("Container", "Population", "Fixed"), 
#'                    data=BeetlesFemale, nboot=3, npermut=3, adjusted=FALSE)
#'                ' 
#' # variance estimation of random effects, residual and overdispersion 
#' rptPoisson(formula = Egg ~ Treatment + (1|Container) + (1|Population) , 
#'                    grname=c("Container","Population","Residual", "Overdispersion"), 
#'                    data = BeetlesFemale, nboot=3, npermut=3, ratio = FALSE)
#' 
#' @export
#' 

rptPoisson <- function(formula, grname, data, link = c("log", "sqrt"), CI = 0.95, nboot = 1000, 
        npermut = 0, parallel = FALSE, ncores = NULL, ratio = TRUE, adjusted = TRUE) {
        
        # missing values
        no_NA_vals <- stats::complete.cases(data[all.vars(formula)])
        if (sum(!no_NA_vals ) > 0 ){
                warning(paste0(sum(!no_NA_vals), " rows containing missing values were removed"))
                data <- data[no_NA_vals, ]
        } 
        
        # check whether grnames just contain "Residual" or "Overdispersion"
        if (!any((grname != "Residual") & (grname != "Overdispersion") & (grname != "Fixed"))) stop("Specify at least one grouping factor in grname")
        
        # link
        if (length(link) > 1) link <- "log" 
        if (!(link %in% c("log", "sqrt"))) stop("Link function has to be 'log' or 'sqrt'")
        
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
        
        # point estimates of R
        R_pe <- function(formula, data, grname, peYN = FALSE) {
                
                # mod <- suppressWarnings(lme4::glmer(formula = formula, data = data, family = stats::poisson(link = link)))
                mod <- lme4::glmer(formula = formula, data = data, family = stats::poisson(link = link))
                
                # random effect variance data.frame
                VarComps <- as.data.frame(lme4::VarCorr(mod))
                rownames(VarComps) = VarComps$grp
                
                # Check whether Residual is selected
                if (any(grname == "Residual")){
                        output_resid <- TRUE
                        # # delete Residual element
                        grname <- grname[-which(grname == "Residual")]
                }

                # Check whether Fixed is selected
                if (any(grname == "Fixed")){
                        output_fixed <- TRUE
                        # # delete fixed element
                        grname <- grname[-which(grname == "Fixed")]
                }
                
                # groups random effect variances
                var_a <- VarComps[grname, "vcov"]
                names(var_a) <- grname
                
                # intercept on link scale
                beta0 <- unname(lme4::fixef(mod)[1])
                
                # Overdispersion variance
                if (link == "sqrt") {
                        var_r <- VarComps["Overdispersion", "vcov"] + 0.25
                }
                if (link == "log") {
                        estdv <- log(1/exp(beta0)+1)
                        var_r <- VarComps["Overdispersion", "vcov"] + estdv
                }
                
                # Fixed effect variance
                var_f <- stats::var(stats::predict(mod, re.form=NA))

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
                                var_p_link <- sum(VarComps[,"vcov"]) + 0.25
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
                                estdv = log(1/exp(beta0)+1)
                                var_p_link <- sum(VarComps[,"vcov"]) +  estdv
                                if(!adjusted) var_p_link <- var_p_link + var_f
                                R_link = var_a / var_p_link
                                R_r <- var_r / var_p_link
                                R_f_link <- var_f / var_p_link
                                # origial scale
                                EY <- exp(beta0 + (sum(VarComps[,"vcov"]))/2)
                                var_p_org <- EY * (exp(sum(VarComps[,"vcov"])) - 1) + 1
                                if(!adjusted) var_p_org <- var_p_org + var_f
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
        if (nboot > 0)  Ysim <- as.matrix(stats::simulate(mod, nsim = nboot))
       
        bootstr <- function(y, mod, formula, data, grname) {
                data[, names(stats::model.frame(mod))[1]] <- as.vector(y)
                R_pe(formula, data, grname)
        }
        
        warnings_boot <- withWarnings({
                
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
                R_boot <- NA
        }
        })
        # transform bootstrapping repeatabilities into vectors
        boot_org <- as.list(rep(NA, length(grname)))
        boot_link <- as.list(rep(NA, length(grname)))
        if (length(R_boot) == 1) {
                # creating tables when R_boot = NA for simplicity with subsequent processing
                if (is.na(R_boot)) {
                        for(i in c("se_org", "se_link")){
                                assign(i, structure(data.frame(matrix(NA, nrow = length(grname))), 
                                        row.names = grname, names = i), envir = e1)   
                        }
                        for(i in c("CI_org", "CI_link")){
                                assign(i, structure(data.frame(matrix(NA, nrow = length(grname), 
                                        ncol = 2)), row.names = grname), envir = e1)   
                        }
                }
        } else  {
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

        # delete from grname
        if (any(grname == "Residual")){
                output_resid <- TRUE
                # # delete Residual element
                grname <- grname[-which(grname == "Residual")]
        }
        if (any(grname == "Fixed")){
                output_fixed <- TRUE
                # # delete Residual element
                grname <- grname[-which(grname == "Fixed")]
        }

        # significance test by permutation of residuals
        P_permut <- rep(NA, length(grname))
        
        # no permutation test
        if (npermut == 1) {
                R_permut <- NA # earlier: R
                P_permut <- NA
        }
        
        # significance test by permutation of residuals
        # nperm argument just used for parallisation
        
        permut <- function(nperm, formula, mod_red, dep_var, grname, data) {
                if (link == "sqrt") {
                        y_perm <- stats::rpois(nrow(data), 
                        (sqrt(stats::fitted(mod_red)) + sample(stats::resid(mod_red)))^2)
                        }
                if (link == "log") {
                        y_perm <- stats::rpois(nrow(data), 
                        exp(log(stats::fitted(mod_red)) + sample(stats::resid(mod_red))))
                }
                data_perm <- data
                data_perm[dep_var] <- y_perm
                out <- R_pe(formula, data_perm, grname)
                out
        }

        # response variable
        dep_var <- as.character(formula)[2]

        # R_permut <- matrix(rep(NA, length(grname) * npermut), nrow = length(grname))
        P_permut <- structure(data.frame(matrix(NA, nrow = 2, ncol = length(grname)),
                row.names = c("P_permut_org", "P_permut_link")), names = grname)
        
        # for likelihood ratio and permutation test
        terms <- attr(terms(formula), "term.labels")
        randterms <- terms[which(regexpr(" | ", terms, perl = TRUE) > 0)]
        
        warnings_permut <- withWarnings({
                
         if (npermut > 1){
                 for (i in 1:length(grname)) {
                         if (length(randterms) > 1) {
                                 formula_red <- stats::update(formula, eval(paste(". ~ . ", paste("- (1 | ", grname[i], 
                                         ")"))))
                                 mod_red <- lme4::glmer(formula_red, data = data, family = stats::poisson(link = link))
                         } else if (length(randterms) == 1) {
                                 formula_red <- stats::update(formula, eval(paste(". ~ . ", paste("- (", randterms, ")"))))
                                 mod_red <- stats::glm(formula_red, data = data, family = stats::poisson(link = link))
                         }
                 if(parallel == TRUE) {
                         if (is.null(ncores)) {
                                 ncores <- parallel::detectCores()
                                 warning("No core number specified: detectCores() is used to detect the number of \n cores on the local machine")
                         }
                         # start cluster
                         cl <- parallel::makeCluster(ncores)
                         parallel::clusterExport(cl, "R_pe", envir=environment())
                         R_permut <- parallel::parLapply(cl, 1:(npermut-1), permut, formula=formula, 
                                 mod_red=mod_red, dep_var=dep_var, grname=grname, data = data)
                         parallel::stopCluster(cl)
                         
                 } else if (parallel == FALSE) {
                         R_permut <- lapply(1:(npermut - 1), permut, formula, mod_red, dep_var, grname, data)
                 }
                 
                 # adding empirical rpt 
                 R_permut <- c(list(R), R_permut)
                 }
         }
        })
        
        # equal to boot
        permut_org <- as.list(rep(NA, length(grname)))
        permut_link <- as.list(rep(NA, length(grname)))
        
        if (!(npermut == 1)){
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
        
        for (i in c("LRT_P", "LRT_D", "LRT_red")) assign(i, rep(NA, length(grname)))
        
        for (i in 1:length(grname)) {
                formula_red <- stats::update(formula, eval(paste(". ~ . ", paste("- (1 | ", grname[i], 
                        ")"))))
                LRT_red[i] <- as.numeric(stats::logLik(lme4::glmer(formula = formula_red, data = data, 
                        family = stats::poisson(link = link))))
                LRT_D[i] <- as.numeric(-2 * (LRT_red[i] - LRT_mod))
                LRT_P[i] <- ifelse(LRT_D[i] <= 0, 1, stats::pchisq(LRT_D[i], 1, lower.tail = FALSE)/2)
                # LR <- as.numeric(-2*(logLik(lme4::lmer(update(formula, eval(paste('. ~ . ',
                # paste('- (1 | ', grname[i], ')') ))), data=data))-logLik(mod))) P.LRT[i] <-
                # ifelse(LR<=0, 1, pchisq(LR,1,lower.tail=FALSE)/2)
        }
  
        P <- cbind(LRT_P, t(P_permut))
        row.names(P) <- grname
        
        # add Residual = NA for S3 functions to work
        if(any(grname_org == "Residual")){
                # grname <- grname_org
                P <- rbind(P, NA)
                row.names(P)[nrow(P)] <- "Residual"
                permut_link$Residual <- rep(NA, length(permut_link[[1]]))
                permut_org$Residual <- rep(NA, length(permut_org[[1]]))
        }
        # add Fixed = NA for S3 functions to work
        if(any(grname_org == "Fixed")){
                # grname <- grname_org
                P <- rbind(P, NA)
                row.names(P)[nrow(P)] <- "Fixed"
                permut_link$Fixed <- rep(NA, length(permut_link[[1]]))
                permut_org$Fixed <- rep(NA, length(permut_org[[1]]))
        }
        
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
                LRT = list(LRT_mod = LRT_mod, LRT_red = LRT_red, LRT_D = LRT_D, LRT_df = LRT_df, 
                LRT_P = LRT_P), 
                ngroups = unlist(lapply(data[grname], function(x) length(unique(x)))), 
                nobs = nrow(data), mod = mod, ratio = ratio, 
                all_warnings = list(warnings_boot = warnings_boot, warnings_permut = warnings_permut))
        
        class(res) <- "rpt"
        return(res)
} 
