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



#' Bootstrapping for non-gaussian functions
#' 
#' @param bootstr bootstrap function. Re-assigns response simulated by simulate.merMod to data and
#' estimates R with the R_pe function.
#' @param R_pe Function to estimate Repeatabilities and Variances for grouping factors,
#' Residuals, Overdispersion and Fixed effects. 
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

bootstrap_nongaussian <- function(bootstr, R_pe, formula, data, Ysim, mod, grname, grname_org, nboot, parallel, ncores, CI) {
        
        e_boot <- environment()
        
        warnings_boot <<- with_warnings({
                
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
                        R_boot <- unname(lapply(Ysim, bootstr, mod, formula, data , 
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
                        boot_org[[i]] <- unlist(lapply(R_boot, function(x) x["R_org", grname_org[i]]))
                        boot_link[[i]] <- unlist(lapply(R_boot, function(x) x["R_link", grname_org[i]]))
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


