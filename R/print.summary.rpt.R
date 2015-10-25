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
#' # repeatability estimation for weight (body mass)
#' data(BodySize)
#' rpt.Weight <- rpt.mcmcLMM(data = BodySize, Weight, BirdID) 
#' summary(rpt.Weight)
#' 
print.summary.rpt <- function(x, ...) {
    
    if (x$datatype == "Gaussian" & x$method == "corr") {
        cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n", 
            "Call = ", deparse(x$call), "\n", "Data: ", x$nobs, " observations in ", 
            x$ngroups, " groups", "\n\n", sep = "")
        cat("Repeatability:", "\n")
        print(format(x$rpt, digits = 3), row.names = FALSE)
        cat("\n")
        cat("Bootstrapping and Permutation test:", "\n")
        print(format(rbind(x$boot, x$permut), digits = 3))
    }
    
    if (x$datatype == "Gaussian" & x$method == "ANOVA") {
        cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n", 
            "Call = ", deparse(x$call), "\n", "Data: ", x$nobs, " observations in ", 
            x$ngroups, " groups", "\n\n", sep = "")
        cat("Repeatability:", "\n")
        print(format(x$rpt, digits = 3, width = 6), row.names = FALSE)
        cat("\n")
        cat("Permutation test:", "\n")
        print(format(x$permut, digits = 3, width = 6), row.names = FALSE)
        cat("\n\n")
        cat("Analysis of Variance Table", "\n")
        cat("Response:", as.character(x$call)[2], "\n")
        print(format(x$anovatab, digits = 3, width = 6))
    }
    
    if (x$datatype == "Gaussian" & x$method == "LMM.REML" & length(x$R) == 1) {
        cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n", 
            "Call = ", deparse(x$call), "\n", "Data: ", x$nobs, " observations in ", 
            x$ngroups, " groups", "\n\n", sep = "")
        cat("Repeatability:", "\n")
        print(format(x$rpt, digits = 3, width = 6), row.names = FALSE)
        cat("\n")
        cat("Bootstrapping and Permutation test:", "\n")
        print(format(rbind(x$boot, x$permut), digits = 3, width = 6))
        cat("\n", "Likelihood ratio test: ", "\n", "logLik full model = ", x$LRT["LRT.mod"], 
            "\n", "logLik red. model = ", x$LRT["LRT.red"], "\n", "D  = ", signif(x$LRT["LRT.D"], 
                3), ", ", "df = ", x$LRT["LRT.df"], ", ", "P  = ", signif(x$LRT["LRT.P"], 
                3), sep = "")
    }
    
    if (x$datatype == "Gaussian" & length(x$P) > 1 & length(x$R) > 1) {
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
                ", ", "P  = ", signif((x$LRT[["LRT.P"]][i]), 3), sep = "")
            cat("\n")
            cat("----------------------------------------")
        }
    }
    
    
    if (x$datatype == "Gaussian" & x$method == "LMM.MCMC") {
        cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n", 
            "Call = ", deparse(x$call), "\n", "Data: ", x$nobs, " observations in ", 
            x$ngroups, " groups", "\n\n", sep = "")
        cat("Repeatability:", "\n")
        print(format(x$rpt, digits = 3), row.names = FALSE)
        cat("\n")
        cat("MCMC parameters:")
        cat("\n")
        cat("Length = ", length(x$R.post), ", ", "Burnin = ", x$MCMCpars[1], ", ", "End = ", 
            x$MCMCpars[2], ", ", "Thinning = ", x$MCMCpars[3], sep = "")
    }
    
    if (x$datatype != "Gaussian" & x$method == "PQL") {
        cat("\n", "Repeatability calculation using the ", x$method, " method and ", x$link, 
            "link", "\n\n", "Call = ", deparse(x$call), "\n", "Data: ", x$nobs, " observations in ", 
            x$ngroups, " groups", "\n\n", "Estimated overdispersion (omega) = ", x$omega, 
            "\n\n", sep = "")
        cat("Repeatability:", "\n")
        print(format(x$rpt, digits = 3))
        cat("\n")
        cat("Bootstrapping:", "\n")
        print(format(x$boot, digits = 3))
        cat("\n")
        cat("Permutation:", "\n")
        print(format(x$perm, digits = 3))
    }
    
    if (x$datatype != "Gaussian" & x$method == "MCMC") {
        # x$datatype!='Gaussian' &
        cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n", 
            "Call = ", deparse(x$call), "\n", "Data: ", x$nobs, " observations in ", 
            x$ngroups, " groups", "\n\n", sep = "")
        cat("Repeatability:", "\n")
        print(format(x$rpt, digits = 3))
        cat("\n")
        cat("MCMC parameters:")
        cat("\n")
        cat("Length = ", length(x$R.post$R.link), ", ", "Burnin = ", x$MCMCpars[1], ", ", 
            "End = ", x$MCMCpars[2], ", ", "Thinning = ", x$MCMCpars[3], sep = "")
        # print(format(data.frame('Length' = length(x$R.post$R.link), 'Burnin' =
        # x$MCMCpars[1], 'End' = x$MCMCpars[2],'Thin' = x$MCMCpars[3]), digits = 3, width =
        # 6), row.names = '')
    }
    
} 
