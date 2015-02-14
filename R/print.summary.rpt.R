#' Prints the summary of a rpt object
#'
#' Displays the summary of an rpt object (i.e. the result of a rpt function call) an a nice form.
#' 
#' @param x An rpt object returned from one of the rpt functions
#' @param \dots Additional arguments; none are used in this method.
#'
#' @references 
#' Nakagawa, S. and Schielzeth, H. (2010) \emph{Repeatability for Gaussian and 
#'              non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
#' 
#' @author Holger Schielzeth  (holger.schielzeth@@ebc.uu.se) & 
#'      Shinichi Nakagawa (shinichi.nakagawa@@otago.ac.nz)
#'      
#' @seealso \link{rpt}, \link{rpt.corr}, \link{rpt.aov}, \link{rpt.remlLMM}, \link{rpt.mcmcLMM},
#'          \link{rpt.binomGLMM.add}, \link{rpt.binomGLMM.multi}, \link{rpt.poisGLMM.add}, \link{rpt.poisGLMM.multi}
#' 
#' @examples  
#' # repeatability estimation for weight (body mass)
#' data(BodySize)
#' attach(BodySize)
#' print(rpt.Weight <- rpt.mcmcLMM(Weight, BirdID))
#' detach(BodySize)
#'       
#' @keywords models
#' 
#' @export
#' 
#' 
#' @examples  
#' repeatability estimation for weight (body mass)
#' data(BodySize)
#' attach(BodySize)
#' rpt.Weight <- rpt.mcmcLMM(Weight, BirdID)
#' print(rpt.Weight)  # alternative call to print.rpt() through print()
#' detach(BodySize)
#' 
print.summary.rpt <- function(x) {
        
        if(x$datatype=="Gaussian" & x$method == "corr") {
                cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n",
                    "Call = ", deparse(x$call), 
                    "\n", "Data: ", x$nobs, " observations in ", x$ngroups, " groups", "\n\n",
                    sep = "")  
                cat("\n")
                cat("Repeatability:", "\n")
                print(format(x$rpt, digits = 3), row.names = FALSE)
                cat("\n")
                cat("Bootstrapping and Permutation test:", "\n")
                print(format(rbind(x$boot, x$permut), digits = 3))  
        } 
        
        if(x$datatype=="Gaussian" & x$method == "ANOVA") {
                cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n",
                    "Call = ", deparse(x$call), 
                    "\n", "Data: ", x$nobs, " observations in ", x$ngroups, " groups", "\n\n",
                    sep = "")   
                cat("\n")
                cat("Repeatability:", "\n")
                print(format(x$rpt, digits = 3), row.names = FALSE)
                cat("\n")
                cat("Permutation test:", "\n")
                print(format(x$permut, digits = 3), row.names = FALSE)  
        } 
        
        if(x$datatype=="Gaussian" & x$method == "LMM.REML") { 
                cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n",
                    "Call = ", deparse(x$call), 
                    "\n", "Data: ", x$nobs, " observations in ", x$ngroups, " groups", "\n\n",
                    sep = "")   
                cat("\n")
                cat("Repeatability:", "\n")
                print(format(x$rpt, digits = 3), row.names = FALSE)
                cat("\n")
                cat("Bootstrapping and Permutation test:", "\n")
                print(format(rbind(x$boot, x$permut), digits = 3))  
                cat("\n", "Likelihood ratio test: ", "\n",
                    "logLik full model = ", x$LRT["LRT.mod"], "\n",
                    "logLik red. model = ", x$LRT["LRT.red"], "\n",
                    "D  = ", signif(x$LRT["LRT.D"], 3), ", " , "df = ", x$LRT["LRT.df"], ", ",
                    "P  = ", signif(x$LRT["LRT.P"], 3), 
                    sep = "")
        } 
        
        if(x$datatype!="Gaussian" & x$method=="PQL") {
                cat("\n", "Repeatability calculation using the ", x$method, " method and ", x$link, "link", "\n\n",
                    "Call = ", deparse(x$call), 
                    "\n", "Data: ", x$nobs, " observations in ", x$ngroups, " groups", "\n\n", 
                    "Estimated overdispersion (omega) = ", x$omega, "\n\n",
                    sep = "") 
                cat("\n")
                cat("Repeatability:", "\n")
                print(format(x$rpt, digits = 3)) 
                cat("\n")
                cat("Bootstrapping:", "\n")
                print(format(x$boot, digits = 3))
                cat("\n")
                cat("Permutation:", "\n")
                print(format(x$perm, digits = 3))
        }
	
}	