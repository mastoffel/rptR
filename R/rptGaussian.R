#' LMM-based Repeatability Estimation for Gaussian Data
#' 
#' Estimates the repeatability from a general linear mixed-effects models fitted by restricted maximum likelihood (REML).
#' @inheritParams rpt
#' 
#' @return 
#' Returns an object of class \code{rpt} that is a a list with the following elements: 
#' \item{call}{Function call.}
#' \item{datatype}{Response distribution (here: 'Gaussian').}
#' \item{CI}{Coverage of the confidence interval as specified by the \code{CI} argument.}
#' \item{R}{\code{data.frame} with point estimates for repeatabilities for each grouping factor
#'       of interest.}
#' \item{se}{\code{data.frame} with approximate standard errors (\emph{se}) for repeatabilities. 
#'      Rows repsresent grouping factors of interest. Note that the distribution might not be symmetrical, 
#'      in which case the \emph{se} is less informative.}
#' \item{CI_emp}{\code{data.frame} containing the (empirical) confidence intervals for the repeatabilities 
#'      estiamted based parametric bootstrapping. Each row represents a grouping factor of interest.}
#' \item{P}{\code{data.frame} with p-values based on likelihood-ratio tests
#'      (first column) and permutation tests (second column). Each row represents a grouping factor 
#'      of interest.}
#' \item{R_boot}{Vector(s) of parametric bootstrap samples for \emph{R}. Each \code{list}
#'       element respesents a grouping factor.}
#' \item{R_permut}{Vector(s) of permutation samples for \emph{R}. Each \code{list}
#'       element represents a grouping factor.}
#' \item{LRT}{\code{list} with two elements. (1) The likelihood for the full model and a \code{data.frame} 
#'      called \code{LRT_table} for the reduced model(s), which includes columns
#'      for the respective grouping factor(s), the likelihood(s) \emph{logL_red}, likelihood ratio(s)
#'      \emph{LR_D}, p-value(s) \emph{LRT_P} and degrees of freedom \emph{LRT_df}} 
#' \item{ngroups}{Number of groups for each grouping level.}
#' \item{nobs}{Number of observations.}
#' \item{mod}{Fitted model.}
#' \item{ratio}{Boolean. TRUE, if ratios have been estimated, FALSE, if variances have been estimated}
#' \item{adjusted}{Boolean. TRUE, if estimates are adjusted}
#' \item{all_warnings}{\code{list} with two elements. 'warnings_boot' and 'warnings_permut' contain
#' warnings from the lme4 model fitting of bootstrap and permutation samples, respectively.}
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
#' Nakagawa, S. & Schielzeth, H. (2010) \emph{Repeatability for Gaussian and 
#' non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
#' 
#' @author Holger Schielzeth  (holger.schielzeth@@uni-jena.de),
#'         Shinichi Nakagawa (s.nakagawa@unsw.edu.au) &
#'         Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'      
#' @seealso \link{rpt}
#' 
#' @examples  
#' 
#' data(BeetlesBody)
#' 
#' # Note: nboot and npermut are set to 3 for speed reasons. Use larger numbers
#' # for the real analysis.
#' 
#' # one random effect
#' rpt_est <- rptGaussian(BodyL ~ (1|Population), grname="Population", 
#'                    data=BeetlesBody, nboot=3, npermut=3)
#' 
#' # two random effects
#' rptGaussian(BodyL ~ (1|Container) + (1|Population), grname=c("Container", "Population"), 
#'                    data=BeetlesBody, nboot=3, npermut=3)
#'                    
#' # unadjusted repeatabilities with fixed effects and 
#' # estimation of the fixed effect variance
#' rptGaussian(BodyL ~ Sex + Treatment + Habitat + (1|Container) + (1|Population), 
#'                   grname=c("Container", "Population", "Fixed"), 
#'                   data=BeetlesBody, nboot=3, npermut=3, adjusted=FALSE)
#'                   
#'                   
#' # two random effects, estimation of variance (instead repeatability)
#' R_est <- rptGaussian(formula = BodyL ~ (1|Population) + (1|Container), 
#'             grname= c("Population", "Container", "Residual"),
#'             data=BeetlesBody, nboot=3, npermut=3, ratio = FALSE)
#' 
#' @export
#' 

rptGaussian <- function(formula, grname, data, CI = 0.95, nboot = 1000, 
        npermut = 0, parallel = FALSE, ncores = NULL, ratio = TRUE, adjusted = TRUE,
        rptOutput = NULL, update = FALSE) {
        
        # delete rows with missing values
        no_NA_vals <- stats::complete.cases(data[all.vars(formula)])
        if (sum(!no_NA_vals ) > 0 ){
                warning(paste0(sum(!no_NA_vals), " rows containing missing values were removed"))
                data <- data[no_NA_vals, ]
        } 
        
        # check whether grnames just contain "Residual" or "Overdispersion" or "Fixed"
        if (!any((grname != "Residual") & (grname != "Overdispersion") & (grname != "Fixed"))) stop("Specify at least one grouping factor in grname")
        
        # fit model
        mod <- lme4::lmer(formula, data = data)
        
        # check for random slopes
        VarComps <- lme4::VarCorr(mod)
        # check whether matrix occurs in VarComps
        check_rs <- sum(unlist(lapply(VarComps[grname], function(x) sum(dim(x)) > 2)))
        randomslopes <- FALSE
        if (check_rs > 0) randomslopes <- TRUE
        
        # extract variance components
        # VarComps <- as.data.frame(lme4::VarCorr(mod))
        
        # checks for bootstraps and permutations
        if (nboot == 1) {
                warning("nboot has to be greater than 1 to calculate a CI and has been set to 0")
                nboot <- 0
        }
        if (nboot < 0) nboot <- 0
        if (npermut < 1) npermut <- 1
        
        # save the original grname
        grname_org <- grname
        
        # additional grname components
        output_resid <- FALSE
        output_overdisp <- FALSE
        output_fixed <- FALSE
        
        
        # check whether Residual, Overdispersion or Fixed is selected and if so, remove it
        # from grname vector
        
        for (component in c("Residual", "Overdispersion", "Fixed")) {
                if (any(grname == component)){
                        grname <- grname[-which(grname == component)]
                        if (component == "Residual") output_resid <- TRUE
                        if (component == "Overdispersion") output_overdisp <- TRUE
                        if (component == "Fixed") output_fixed <- TRUE
                }
        }
        
        
        
        
        # point estimates of R or var
        R_pe <- function(formula, data, grname, mod = NULL, resp = NULL) {
                
                if (!is.null(mod)) {
                        mod <- lme4::refit(mod, newresp = resp)
                } else {
                        # model
                        mod <- lme4::lmer(formula, data)   
                }
                
                VarComps <- lme4::VarCorr(mod)
                
                # Residual variance
                var_e <- attr(VarComps, "sc")^2
                names(var_e) <- "Residual"
                
                # Overdispersion 
                var_o <- var_e
                names(var_o) <- "Overdispersion"
                
                # Fixed effect variance
                var_f <- stats::var(stats::predict(mod, re.form=NA))
                names(var_f) <- "Fixed"
                
                # group variances
                group_vars <- function(grname, VarComps, mod){
                        # check whether component is a matrix (--> random slopes)
                        if (sum(dim(VarComps[[grname]])) > 2 ){
                                sigma <- VarComps[[grname]] 
                                # design matrix subsetted for the elements of sigma
                                Z <- stats::model.matrix(mod)[, colnames(sigma)]
                                # average variance across covariate
                                var_grname <- sum(rowSums((Z %*% sigma) * Z))/stats::nobs(mod)
                        } else {
                                var_grname <- as.numeric(VarComps[[grname]])
                        }
                        var_grname
                }
                
                var_a <- unlist(lapply(grname, group_vars, VarComps, mod))
                names(var_a) <- grname
                
                # without random slopes
                # group variances
                # var_a <- as.numeric(VarComps[grname])
                # names(var_a) <- grname
                
                
                # without random slopes
                # denominator variance
                # var_p <- sum(as.numeric(VarComps)) + attr(VarComps, "sc")^2
                var_VarComps <- unlist(lapply(names(VarComps), group_vars, VarComps, mod))
                var_p <- sum(as.numeric(var_VarComps)) + attr(VarComps, "sc")^2
                
                if (!adjusted) var_p <- var_p + var_f
                
                # return variance instead of repeatability
                if (ratio == FALSE) { 
                        R <- as.data.frame(t(var_a))
                        names(R) <- grname
                        
                        # if residual is selected, add residual variation to the output
                        if (output_resid){
                                R$Residual <- var_e
                        } 
                        if (output_overdisp){
                                R$Overdispersion <- var_o
                        }
                        if (output_fixed){
                                R$Fixed <- var_f
                        }
                        return(R)
                }
                
                # return repeatability
                if (ratio == TRUE) { 
                        
                        R <- var_a/var_p
                        R <- as.data.frame(t(R))
                        names(R) <- grname
                        
                        # check whether to give out Residual
                        if(output_resid){
                                R$Residual <- var_e / var_p
                        }
                        # check whether to give out Overdispersion
                        if(output_overdisp){
                                R$Overdispersion <- var_e / var_p
                        }
                        # check whether to give out Fixed
                        if(output_fixed){
                                R$Fixed <- var_f / var_p
                        }
                        return(R)
                }
        }
        
        R <- R_pe(formula, data, grname) # no bootstrap skipping at the moment
        
        # confidence interval estimation by parametric bootstrapping
        # simulate matrix from which to draw y
        if (nboot > 0)  Ysim <- as.matrix(stats::simulate(mod, nsim = nboot))
        
        # bootstrapping function
        bootstr <- function(y, mod, formula, data, grname) {
                # data[, names(stats::model.frame(mod))[1]] <- as.vector(y)
                resp <- as.vector(y)
                # add mod and resp to use lme4::refit instead of fitting a new model
                R_pe(formula, data, grname, mod = mod, resp = resp)
        }
        
        num_iter <- NULL
        
        warnings_boot <- with_warnings({
                
                if (nboot > 0 & parallel == TRUE) {
                        if (is.null(ncores)) {
                                ncores <- parallel::detectCores() - 1
                                warning("No core number specified: detectCores() is used to detect the number of \n cores on the local machine")
                        }
                        # start cluster
                        cl <- parallel::makeCluster(ncores)
                        parallel::clusterExport(cl, "R_pe", envir=environment())
                        R_boot <- unname(parallel::parApply(cl = cl, Ysim, 2, bootstr, mod = mod, formula = formula, 
                                data = data, grname = grname))
                        parallel::stopCluster(cl)
                }
                
                if (nboot > 0 & parallel == FALSE) {
                        
                        cat("Bootstrap Progress:\n")
                        R_boot <- unname(pbapply::pbapply(Ysim, 2, bootstr, mod = mod, formula = formula, data = data, 
                                grname = grname))
                        
                }
                if (nboot == 0) {
                        R_boot <- NA
                }
                
        })
        
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
                boot <- do.call(rbind, R_boot)
                
                ## addition
                if (update){
                        if (is.null(rptOutput)) stop("provide rpt object for rptOutput argument")
                        boot <- rbind(boot, rptOutput$R_boot)
                }
                
                CI_emp <- as.data.frame(t(apply(boot, 2, calc_CI)))
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
        permut <- function(nperm, formula, data, mod_red, dep_var, grname, i, mod) {
                y_perm <- stats::fitted(mod_red) + sample(stats::resid(mod_red))
                data_perm <- data
                data_perm[dep_var] <- y_perm
                
                # add mod and resp to use lme4::refit instead of fitting a new model
                out <- R_pe(formula, data_perm, grname[i], mod = mod, resp = y_perm)
                out
        }
        
        dep_var <- as.character(formula)[2]
        # one random effect, uses stats::lm()
        # multiple random effects, uses lmer()
        
        R_permut <- data.frame(matrix(rep(NA, length(grname) * npermut), nrow = length(grname)))
        P_permut <- rep(NA, length(grname))
        
        if (update){
                if (is.null(rptOutput)) stop("provide rpt object for rptOutput argument")
                
                old_perm_length <- length(rptOutput$R_permut[[1]])
                R_permut <- data.frame(matrix(rep(NA, length(grname) * (npermut + old_perm_length)), nrow = length(grname)))
                # add new R permut without empirical point estimate to old R permut
                if (npermut > 0) npermut <- npermut + 1 # account for deleting the point estimate in the update
        }
        
        # function for the reduced model in permut and LRT tests
        mod_fun <- ifelse(length(randterms) == 1, stats::lm, lme4::lmer)
        
        warnings_permut <- with_warnings({
                
                if (npermut > 1){
                        for (i in 1:length(grname)) {
                                # from random terms
                                randterm <-  randterms[grep(grname[i], randterms)]
                                # formula_red <- stats::update(formula, eval(paste(". ~ . ", paste("- (1 | ", grname[i], ")")))) ## check that random slopes work
                                formula_red <- stats::update(formula, eval(paste(". ~ . ", paste("- (", randterm, ")")))) ## check that random slopes work
                                mod_red <- mod_fun(formula_red, data = data)
                                if(parallel == TRUE) {
                                        if (is.null(ncores)) {
                                                ncores <- parallel::detectCores()
                                                warning("No core number specified: detectCores() is used to detect the number of \n cores on the local machine")
                                        }
                                        # start cluster
                                        cl <- parallel::makeCluster(ncores)
                                        parallel::clusterExport(cl, "R_pe", envir=environment())
                                        R_permut[i, ] <- c(R[i], as.numeric(unlist(parallel::parSapply(cl, 1:(npermut-1), permut, formula, data, mod_red, dep_var, grname, i, mod))))
                                        parallel::stopCluster(cl)
                                } else if (parallel == FALSE) {
                                        cat("Permutation Progress for", grname[i], ":\n")
                                        R_permut[i, ] <- c(R[i], as.numeric(unlist(pbapply::pbreplicate(npermut - 1, permut(formula=formula, data = data, 
                                                mod_red=mod_red, dep_var=dep_var, grname=grname, i=i, mod = mod)))))
                                        
                                }
                                
                                ## addition
                                if (update){
                                        if (is.null(rptOutput)) stop("provide rpt object for rptOutput argument")
                                        
                                        R_permut[i, ] <- c(rptOutput$R_permut[[i]], unlist(R_permut[i, ])[-1])
                                        # add new R permut without empirical point estimate to old R permut
                                }
                                
                                P_permut[i] <- sum(R_permut[i, ] >= unlist(R[i]))/npermut
                        }
                }
                
        })
        # name R_permut and P_permut
        row.names(R_permut) <- grname
        names(P_permut) <- grname
        
        
        ## likelihood-ratio-test
        LRT_mod <- as.numeric(stats::logLik(mod))
        LRT_df <- rep(1, length(grname)) # (2*2) dimension is 3 df, 3*3 is 6 df , if k is dimension of matrix:
                                         # 3 free parameters random intercept, random slop and correlation
        
        # preassign
        for (i in c("LRT_P", "LRT_D", "LRT_red")) assign(i, rep(NA, length(grname)))
        # function
        # mod_fun <- ifelse(length(randterms) == 1, stats::lm, lme4::lmer)
        
        for (i in 1:length(grname)) {
                randterm <-  randterms[grep(grname[i], randterms)]
                # formula_red <- stats::update(formula, eval(paste(". ~ . ", paste("- (1 | ", grname[i], ")")))) ## check that random slopes work
                formula_red <- stats::update(formula, eval(paste(". ~ . ", paste("- (", randterm, ")")))) ## check that random slopes work
                LRT_red[i] <- as.numeric(stats::logLik(mod_fun(formula = formula_red, data = data)))
                LRT_D[i] <- as.numeric(-2 * (LRT_red[i] - LRT_mod))
                LRT_P[i] <- ifelse(LRT_D[i] <= 0, 1, stats::pchisq(LRT_D[i], 1, lower.tail = FALSE)/2) # instead of 1: LRT_df[i]
        }
        LRT_table <- data.frame(logL_red = LRT_red, LR_D = LRT_D, LRT_P = LRT_P, LRT_df =  LRT_df, stringsAsFactors = FALSE)
        row.names(LRT_table) <- grname
        
        P <- as.data.frame(cbind(LRT_P, P_permut))
        row.names(P) <- grname
        
        for (component in c("Residual", "Overdispersion", "Fixed")) {
                if(any(grname_org == component)){
                        P <- rbind(P, as.numeric(NA))
                        row.names(P)[nrow(P)] <- component
                        R_permut <- rbind(R_permut, as.numeric(NA))
                        row.names(R_permut)[nrow(R_permut)] <- component
                        LRT_table <- rbind(LRT_table, as.numeric(NA))
                        row.names(LRT_table)[nrow(LRT_table)] <- component
                }
        }
        
        
        res <- list(call = match.call(), 
                datatype = "Gaussian", 
                CI = CI, 
                R = R, 
                se = se,
                CI_emp = CI_emp, 
                P = P,
                R_boot = boot, 
                R_permut = lapply(as.data.frame(t(R_permut)), function(x) return(x)),
                LRT = list(LRT_mod = LRT_mod, LRT_table = LRT_table), 
                ngroups = unlist(lapply(data[grname], function(x) length(unique(x)))), 
                nobs = nrow(data), mod = mod, ratio = ratio, adjusted = adjusted,
                all_warnings = list(warnings_boot = warnings_boot, warnings_permut = warnings_permut))
        
        class(res) <- "rpt"
        return(res)
} 
