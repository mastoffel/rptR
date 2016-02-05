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
#' @author Holger Schielzeth  (holger.schielzeth@@ebc.uu.se),
#'         Shinichi Nakagawa (shinichi.nakagawa@@otago.ac.nz) &
#'         Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'      
#' @seealso \link{rpt}, \link{rpt.corr}, \link{rpt.aov}, \link{rpt.remlLMM}, \link{rpt.mcmcLMM},
#'          \link{rpt.binomGLMM.add}, \link{rpt.binomGLMM.multi}, \link{rpt.poisGLMM.add}, \link{rpt.poisGLMM.multi}
#' 
#' @keywords models
#' 
#' @export
#' 
#' 
#' @examples  
#' 
print.summary.rpt <- function(x, ...) {
    
    if (x$datatype == "Poisson" | x$datatype == "Binary" | x$datatype == "Proportion" ) {
            cat("\n", "Repeatability calculation using glmer ",  "\n\n", 
                    "Call = ", deparse(x$call), "\n", "Data: ", x$nobs, " observations", sep = "")
            cat("\n")
            cat("----------------------------------------")
            for (i in 1:ncol(x$R)) {
                    cat("\n\n")
                    cat(names(x$R)[i], " (", x$ngroups[i], " groups)", "\n\n", sep = "")
                    cat("Repeatability overview:", "\n")
                    print(format(x$rpt[[i]], digits = 3, width = 6))
                    cat("\n\n")
                    cat("Bootstrapping:", "\n")
                    print(format(x$boot[[i]], digits = 3, width = 6))
                    cat("\n")
                    cat("Permutation test:", "\n")
                    print(format(x$permut[[i]], digits = 3, width = 6))
                    cat("\n", "Likelihood ratio test: ", "\n", "logLik full model = ", (x$LRT[["LRT_mod"]]), 
                            "\n", "logLik red. model = ", (x$LRT[["LRT_red"]][i]), "\n", "D  = ", 
                            signif((x$LRT[["LRT_D"]][i]), 3), ", ", "df = ", unname(x$LRT[["LRT_df"]]), 
                            ", ", "P_val  = ", signif((x$LRT[["LRT_P"]][i]), 3), sep = "")
                    cat("\n")
                    cat("----------------------------------------")
            }
    }
    
    if (x$datatype == "Gaussian") {
        cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n", 
            "Call = ", deparse(x$call), "\n", "Data: ", x$nobs, " observations", sep = "")
        cat("\n")
        cat("----------------------------------------")
        for (i in 1:length(x$R)) {
            cat("\n\n")
            cat(names(x$R)[i], " (", x$ngroups[i], " groups)", "\n\n", sep = "")
            cat("Repeatability:", "\n")
            print(format(x$rpt[[i]], digits = 3, width = 6), row.names = FALSE)
            cat("\n")
            cat("Bootstrapping and Permutation test:", "\n")
            print(format(rbind(x$boot[[i]], x$permut[[i]]), digits = 3, width = 6))
            cat("\n", "Likelihood ratio test: ", "\n", "logLik full model = ", (x$LRT[["LRT.mod"]]), 
                "\n", "logLik red. model = ", (x$LRT[["LRT.red"]][i]), "\n", "D  = ", 
                signif((x$LRT[["LRT.D"]][i]), 3), ", ", "df = ", unname(x$LRT[["LRT.df"]]), 
                ", ", "P_val  = ", signif((x$LRT[["LRT.P"]][i]), 3), sep = "")
            cat("\n")
            cat("----------------------------------------")
        }
    }
    
} 
