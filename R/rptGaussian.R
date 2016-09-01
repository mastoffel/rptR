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
#' # two random effects, estimation of variance (instead repeatability)
#' rptGaussian(formula = BodyL ~ (1|Population) + (1|Container), 
#'             grname= c("Population", "Container", "Residual"),
#'             data=BeetlesBody, nboot=3, npermut=3, ratio = FALSE)
#' 
#' @export
#' 

rptGaussian <- function(formula, grname, data, CI = 0.95, nboot = 1000, 
        npermut = 0, parallel = FALSE, ncores = NULL, ratio = TRUE, adjusted = TRUE) {
        
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
        R_pe <- function(formula, data, grname) {
                
                # model
                mod <- lme4::lmer(formula, data)
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
                var_a <- as.numeric(VarComps[grname])
                names(var_a) <- grname
                
                # denominator variance
                var_p <- sum(as.numeric(VarComps)) + attr(VarComps, "sc")^2
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
                data[, names(stats::model.frame(mod))[1]] <- as.vector(y)
                R_pe(formula, data, grname)
        }
        
        warnings_boot <- .with_warnings({
                
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
        
        # function for the reduced model in permut and LRT tests
        mod_fun <- ifelse(length(randterms) == 1, stats::lm, lme4::lmer)
        
        warnings_permut <- .with_warnings({
                
        if (npermut > 1){
                for (i in 1:length(grname)) {
                                formula_red <- stats::update(formula, eval(paste(". ~ . ", paste("- (1 | ", grname[i], ")"))))
                                mod_red <- mod_fun(formula_red, data = data)
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
                                        mod_red=mod_red, dep_var=dep_var, grname=grname, i=i)))))
                                P_permut[i] <- sum(R_permut[i, ] >= unlist(R[i]))/npermut
                        }
                }
        }
                
        })
        # name R_permut and P_permut
        row.names(R_permut) <- grname
        names(P_permut) <- grname
  
        
        ## likelihood-ratio-test
        LRT_mod <- as.numeric(stats::logLik(mod))
        LRT_df <- rep(1, length(grname))
        
        # preassign
        for (i in c("LRT_P", "LRT_D", "LRT_red")) assign(i, rep(NA, length(grname)))
        # function
        # mod_fun <- ifelse(length(randterms) == 1, stats::lm, lme4::lmer)
        
        for (i in 1:length(grname)) {
                formula_red <- stats::update(formula, eval(paste(". ~ . ", paste("- (1 | ", grname[i], ")"))))
                LRT_red[i] <- as.numeric(stats::logLik(mod_fun(formula = formula_red, data = data)))
                LRT_D[i] <- as.numeric(-2 * (LRT_red[i] - LRT_mod))
                LRT_P[i] <- ifelse(LRT_D[i] <= 0, 1, stats::pchisq(LRT_D[i], 1, lower.tail = FALSE)/2)
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
