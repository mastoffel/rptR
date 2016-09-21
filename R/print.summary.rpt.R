#' Prints the summary of a rpt object
#'
#' Displays the summary of an rpt object (i.e. the result of a rpt function call) in an 
#' extended form.
#' 
#' @param x An rpt object returned from one of the rpt functions
#' @param \dots Additional arguments; none are used in this method.
#'
#' @references 
#' Nakagawa, S. and Schielzeth, H. (2010) \emph{Repeatability for Gaussian and 
#' non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
#' 
#' @author Holger Schielzeth  (holger.schielzeth@@uni-jena.de), 
#'         Shinichi Nakagawa (s.nakagawa@unsw.edu.au),
#'         Martin Stoffel (martin.adam.stoffel@@gmail.com) 
#'      
#' @keywords models
#' 
#' @export
#' 
#' 
#' 
print.summary.rpt <- function(x, ...) {
        
        if (x$ratio == TRUE) {
                header_gaussian <- "Repeatability estimation using the lmm method"
                header_nongaussian <- "Repeatability estimation using glmer method"
                PE <- "Repeatability estimation overview:"
        } else if (x$ratio == FALSE) {
                header_gaussian <- "Variance estimation using the lmm method"
                header_nongaussian <- "Variance estimation using glmer method"
                PE <- "Variance estimation overview:"
                x$rpt <- lapply(x$rpt, function(x){
                        names(x)[1] <- "Var"
                        x
                })
        }
        
    
    if (x$datatype == "Poisson" | x$datatype == "Binary" | x$datatype == "Proportion" ) {
            cat("\n",  header_nongaussian,  "\n\n", 
                    "Call = ", gsub("^\\s+", "", deparse(x$call)), "\n\n", "Data: ", x$nobs, " observations", sep = "")
            cat("\n")
            cat("----------------------------------------")
            for (i in 1:ncol(x$R)) {
                    cat("\n\n")
                    if (!is.na(x$ngroups[i])) {
                            cat(names(x$R)[i], " (", x$ngroups[i], " groups)", "\n\n", sep = "")
                    } else {
                            cat(names(x$R)[i] , "\n\n",sep = "")
                    }
                    cat(PE, "\n")
                    print(format(x$rpt[[i]], digits = 3, width = 6))
                    cat("\n\n")
                    cat("Bootstrapping:", "\n")
                    print(format(x$boot[[i]], digits = 3, width = 6))
                    cat("\n")
                    cat("Permutation test:", "\n")
                    print(format(x$permut[[i]], digits = 3, width = 6))
                    cat("\n", "Likelihood ratio test: ", "\n", "logLik full model = ", (x$LRT[["LRT_mod"]]), 
                            "\n", "logLik red. model = ", (x$LRT$LRT_table[["logL_red"]][i]), "\n", "D  = ", 
                            signif((x$LRT$LRT_table[["LR_D"]][i]), 3), ", ", "df = ", x$LRT$LRT_table[["LRT_df"]][i] , 
                            ", ", "P = ", signif((x$LRT$LRT_table[["LRT_P"]][i]), 3), sep = "")
                    cat("\n\n")
                    cat("----------------------------------------")
                    cat("\n")
            }
    }
    
    if (x$datatype == "Gaussian") {
        cat("\n", header_gaussian, "\n\n", 
            "Call = ", gsub("^\\s+", "", deparse(x$call)) , "\n\n", "Data: ", x$nobs, " observations", sep = "")
        cat("\n")
        cat("----------------------------------------")
        for (i in 1:length(x$R)) {
            cat("\n\n")
            if (!is.na(x$ngroups[i])) {
                cat(names(x$R)[i], " (", x$ngroups[i], " groups)", "\n\n", sep = "")
            } else {
                cat(names(x$R)[i] , "\n\n",sep = "")
            }
            cat(PE, "\n")
            print(format(x$rpt[[i]], digits = 3, width = 6), row.names = FALSE)
            cat("\n")
            cat("Bootstrapping and Permutation test:", "\n")
            print(format(rbind(x$boot[[i]], x$permut[[i]]), digits = 3, width = 6))
            cat("\n", "Likelihood ratio test: ", "\n", "logLik full model = ", (x$LRT[["LRT_mod"]]), 
                    "\n", "logLik red. model = ", (x$LRT$LRT_table[["logL_red"]][i]), "\n", "D  = ", 
                    signif((x$LRT$LRT_table[["LR_D"]][i]), 3), ", ", "df = ", x$LRT$LRT_table[["LRT_df"]][i], 
                    ", ", "P = ", signif(x$LRT$LRT_table[["LRT_P"]][i], 3), sep = "")
            cat("\n\n")
            cat("----------------------------------------")
            cat("\n")
        }
    }
    
} 
