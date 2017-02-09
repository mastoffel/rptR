#' Captures and suppresses (still to find out why) warnings of an expression
#'
#' This function is used within rptR to capture lme4 model fitting warnings in the
#' bootstrap and permutation procedures.
#'
#' @param expr An expression, such as the sequence of code used by rptR to calculate
#' bootstrap or permutation estimates
#' @keywords internal


with_warnings <- function(expr) {
        myWarnings <- NULL
        wHandler <- function(w) {
                myWarnings <<- c(myWarnings, list(w))
                invokeRestart("muffleWarning")
        }
        val <- withCallingHandlers(expr, warning = wHandler)
        list(warnings = myWarnings)
} 

#' Calculates / extracts variance components from random effects and random slopes
#'
#' This function uses the method from Paul Johnson to compute the average group
#' variance across the levels of a covariate.
#'
#' @param grname The name of a grouping factor, usually accessed by looping over the
#' grname argument of the rptR functions.
#' @param VarComps A list. Output of the lme4::VarCorr function.
#' @param mod An lme4 model object.
#' @keywords internal
#' 
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


#' Bootstrapping for non-gaussian functions (internal use)
#' 
#' @param bootstr bootstrap function. Re-assigns response simulated by simulate.merMod to data and estimates R with the R_pe function.
#' @param R_pe Function to estimate Repeatabilities and Variances for grouping factors, Residuals, Overdispersion and Fixed effects. 
#' @param formula lme4 model formula
#' @param data data.frame given as original input
#' @param Ysim data.frame with simulated response variables from simulate.merMod
#' @param mod fitted lme4 model
#' @param grname original grnames vector without Residual or Fixed
#' @param grname_org original grnames vector
#' @param nboot number of bootstraps, equal to columns in Ysim
#' @param parallel boolean
#' @param ncores number of cores specified, defaults to NULL
#' @param CI confidence interval, defaults to 0.95
#' @keywords internal

bootstrap_nongaussian <- function(bootstr, R_pe, formula, data, Ysim, mod, grname, grname_org, nboot, parallel, ncores, CI, rptObj, update) {
        
        e_boot <- environment()
        
        warnings_boot <- with_warnings({
                
                # to do: preallocate R_boot
                if (nboot > 0 & parallel == TRUE) {
                        if (is.null(ncores)) {
                                ncores <- parallel::detectCores() - 1
                                warning("No core number specified: detectCores() is used to detect the number of \n cores on the local machine")
                        }
                        # start cluster
                        cl <- parallel::makeCluster(ncores)
                        parallel::clusterExport(cl, "R_pe", envir=environment())
                        R_boot <- unname(parallel::parLapply(cl, Ysim, bootstr, mod, formula, 
                                data, grname))
                        parallel::stopCluster(cl)
                }
                if (nboot > 0 & parallel == FALSE) {
                        cat("Bootstrap Progress:\n")
                        R_boot <- unname(pbapply::pblapply(Ysim, bootstr, mod, formula, data , 
                                grname))
                }
                if (nboot == 0) {
                        R_boot <- NA
                }
                
        })
        
        # transform bootstrapping repeatabilities into vectors
        boot_org <- as.list(rep(NA, length(grname_org)))
        boot_link <- as.list(rep(NA, length(grname_org)))
        if (length(R_boot) == 1) {
                # creating tables when R_boot = NA
                if (is.na(R_boot)) {
                        # for(i in c("CI_org", "CI_link", "se_org", "se_link")) assign(i, NA, envir = e1)
                        for(i in c("se_org", "se_link")){
                                assign(i, structure(data.frame(matrix(NA, 
                                        nrow = length(grname_org))), row.names = grname_org, names = i), envir = e_boot)   
                        }
                        for(i in c("CI_org", "CI_link")){
                                assign(i, structure(data.frame(matrix(NA, nrow = length(grname_org), ncol = 2)), 
                                        row.names = grname_org), envir = e_boot)   
                        }
                        
                }
        } else {
                for (i in 1:length(grname_org)) {
                        
                        if (update){
                                if (is.null(rptObj)) stop("provide rpt object for rptObj argument")
                                boot_org[[i]] <- c(rptObj$R_boot_org[[i]], unlist(lapply(R_boot, function(x) x["R_org", grname_org[i]])))
                                boot_link[[i]] <- c(rptObj$R_boot_link[[i]],unlist(lapply(R_boot, function(x) x["R_link", grname_org[i]])))
                        } else {
                                boot_org[[i]] <- unlist(lapply(R_boot, function(x) x["R_org", grname_org[i]]))
                                boot_link[[i]] <- unlist(lapply(R_boot, function(x) x["R_link", grname_org[i]])) 
                        }
                        
         
                }
                names(boot_org) <- grname_org
                names(boot_link) <- grname_org
                
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

        
        out <- list(R_boot = R_boot, boot_org = boot_org, boot_link = boot_link, CI_org = CI_org, CI_link = CI_link,
                se_org = se_org, se_link = se_link)

}




#' Permutation function for non-gaussian functions (internal use)
#' 
#' @inheritParams bootstrap_nongaussian
#' @param permut permutation function which permutes residuals and calculates R
#' @param dep_var original response variable
#' @param link link function
#' @param family respnse family (so far just binomial or poisson)
#' @param npermut number of permutations
#' @param R point estimate to concetenate with permutations
#' 
#' @keywords internal

permut_nongaussian <- function(permut, R_pe, formula, data, dep_var, grname, npermut, parallel, 
                               ncores, link, family, R, rptObj, update){
        
  
        # predefine if no permutation test
        if (npermut == 1) {
                R_permut <- NA
                P_permut <- NA
        }
        
        # R_permut <- matrix(rep(NA, length(grname) * npermut), nrow = length(grname))
        P_permut <- structure(data.frame(matrix(NA, nrow = 2, ncol = length(grname)),
                row.names = c("P_permut_org", "P_permut_link")), names = grname)
        
        # for likelihood ratio and permutation test
        terms <- attr(terms(formula), "term.labels")
        randterms <- terms[which(regexpr(" | ", terms, perl = TRUE) > 0)]
        
        # at the moment its not possible to fit the same grouping factor in more than one random effect
        check_modelspecs <- sapply(grname, function(x) sum(grepl(x, randterms)))
        if (any(check_modelspecs > 1)){
                stop("Fitting the same grouping factor in more than one random 
                        effect terms is not possible at the moment")  
        } 
        
        # function for the reduced model in permut and LRT tests
        mod_fun <- ifelse(length(randterms) == 1, stats::glm, lme4::glmer)
        
        # family
        if (family == "poisson") family_fun <- stats::poisson
        if (family == "binomial") family_fun <- stats::binomial
   
        e_permut <- environment()
        
        if (update){
                if (is.null(rptObj)) stop("provide rpt object for rptObj argument")
                # one more permutation as we don't add the empirical point estimate again
                if (npermut > 0)  npermut <- npermut + 1
        }
        
        warnings_permut <- with_warnings({
                
                if (npermut > 1){
                        for (i in 1:length(grname)) {
                                # from random terms
                                randterm <-  randterms[grep(grname[i], randterms)]
                                # formula_red <- stats::update(formula, eval(paste(". ~ . ", paste("- (1 | ", grname[i], ")")))) ## check that random slopes work
                                formula_red <- stats::update(formula, eval(paste(". ~ . ", paste("- (", randterm, ")")))) ## check that random slopes work
                                mod_red <- mod_fun(formula_red, data = data, family = family_fun(link = link))
                                
                                if (!(grname[i] == "Overdispersion")){
                                        
                                if(parallel == TRUE) {
                                        if (is.null(ncores)) {
                                                ncores <- parallel::detectCores()
                                                warning("No core number specified: detectCores() is used to detect the number of \n cores on the local machine")
                                        }
                                        # start cluster
                                        cl <- parallel::makeCluster(ncores)
                                        parallel::clusterExport(cl, "R_pe", envir=environment())
                                        out_permut <- parallel::parLapply(cl, 1:(npermut-1), permut, formula=formula, 
                                                mod=mod_red, dep_var=dep_var, grname=grname[i], data = data)
                                        parallel::stopCluster(cl)
                                        
                                } else if (parallel == FALSE) {
                                        cat(paste0("Permutation Progress for ", grname[i], ":\n"))
                                        out_permut <- pbapply::pblapply(1:(npermut - 1), permut, formula, mod_red, dep_var, grname[i], data)
                                }
                                
                                } 
                                # in case its overdispersion, just add NA data.frames. 
                                if (grname[i] == "Overdispersion") {
                                        out_permut <- lapply(1:(npermut-1), function(x) data.frame("Overdispersion" = c(NA,NA), row.names = c("R_org", "R_link")))
                                }
                                
                                # adding empirical rpt 
                                if(!exists("R_permut", envir = e_permut)) {
                                        R_permut <- out_permut
                                } else {
                                        R_permut <- lapply(1:length(R_permut), function(x) cbind(R_permut[[x]], out_permut[[x]]))
                                }
                                
                        
                        }
                        
                        if (!update){
                                # if we are updated, we don't add the point estimate again
                                R_permut <- c(list(R), R_permut)
                        }
                        
                        # R_permut <- c(list(R), R_permut)
                        R_permut[[1]]["Overdispersion"] <- NA
                }
        })
        #list(R),
        # equal to boot
        permut_org <- as.list(rep(NA, length(grname)))
        permut_link <- as.list(rep(NA, length(grname)))
        
        if (!(length(R_permut) == 1)){
                
                for (i in 1:length(grname)) {
                        
                        if (update){
                                if (is.null(rptObj)) stop("provide rpt object for rptObj argument")
                                permut_org[[i]]<- c(rptObj$R_permut_org[[i]], unlist(lapply(R_permut, function(x) x["R_org", grname[i]])))
                                permut_link[[i]] <- c(rptObj$R_permut_link[[i]], unlist(lapply(R_permut, function(x) x["R_link", grname[i]])))
                        } else {
                                permut_org[[i]] <- unlist(lapply(R_permut, function(x) x["R_org", grname[i]]))
                                permut_link[[i]] <- unlist(lapply(R_permut, function(x) x["R_link", grname[i]]))
                        }
                        
                      
                }
                names(permut_org) <- grname
                names(permut_link) <- grname
        }
        
        
        P_permut["P_permut_org", ] <- unlist(lapply(permut_org, function(x) sum(x >= x[1])))/npermut
        P_permut["P_permut_link", ] <- unlist(lapply(permut_link, function(x) sum(x >= x[1])))/npermut
        names(P_permut) <- names(permut_link)
        
        
        out <- list(P_permut = P_permut, permut_org = permut_org, 
                permut_link = permut_link, warnings_permut = warnings_permut)

}



#' Likelihood ratio test for non-gaussian functions (internal use)
#' 
#' @inheritParams permut_nongaussian
#' 
#' @keywords internal


LRT_nongaussian <- function(formula, data, grname, mod, link, family){
        
        # family
        if (family == "poisson") family_fun <- stats::poisson
        if (family == "binomial") family_fun <- stats::binomial
        
        # for likelihood ratio and permutation test
        terms <- attr(terms(formula), "term.labels")
        randterms <- terms[which(regexpr(" | ", terms, perl = TRUE) > 0)]
        
        # function for the reduced model in permut and LRT tests
        mod_fun <- ifelse(length(randterms) == 1, stats::glm, lme4::glmer)
        
        LRT_mod <- as.numeric(stats::logLik(mod))
        
        # calculate df for random slopes
        VarComps <- lme4::VarCorr(mod)
        mat_dims <- unlist(lapply(VarComps[grname], ncol))
        calc_df <- function(k){
                if (k == 1) df <- 1
                if (k > 1){
                        terms <- attr(terms(formula), "term.labels")
                        current_term <- terms[grep(names(k), terms)]
                        if (grep("0", current_term)){
                                df <- (k*(k-1)/2+k) - 1    
                        } else {
                                df <- k*(k-1)/2+k  
                        }
                } 
                df
        }
        LRT_df <- sapply(mat_dims, calc_df) 

        for (i in c("LRT_P", "LRT_D", "LRT_red")) assign(i, rep(NA, length(grname)))

        for (i in 1:length(grname)) {
                # formula_red <- stats::update(formula, eval(paste(". ~ . ", paste("- (1 | ", grname[i], ")"))))
                randterm <-  randterms[grep(grname[i], randterms)]
                formula_red <- stats::update(formula, eval(paste(". ~ . ", paste("- (", randterm, ")")))) ## check that random slopes work
                
                LRT_red[i] <- as.numeric(stats::logLik(mod_fun(formula = formula_red, data = data,
                        family = family_fun(link = link))))
                LRT_D[i] <- as.numeric(-2 * (LRT_red[i] - LRT_mod))
                LRT_P[i] <- ifelse(LRT_D[i] <= 0, 1, stats::pchisq(LRT_D[i], LRT_df[i], lower.tail = FALSE)) 
        }
        
        # division by 2 if LRT_df = 1
        LRT_P <- LRT_P/ifelse(LRT_df==1,2,1)
        
        LRT_table <- data.frame(logL_red = LRT_red, LR_D = LRT_D, LRT_P = LRT_P, LRT_df =  LRT_df, stringsAsFactors = FALSE)
        row.names(LRT_table) <- grname
        
        out <- list(LRT_mod = LRT_mod, LRT_table = LRT_table)
        
}




