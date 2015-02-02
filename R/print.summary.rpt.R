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
                    "\n", "Data: ", x$nobs, " observations and ", x$ngroups, " groups", "\n\n",
                    sep = "")   
                cat("Repeatability: ", "\n", 
                    "R  = ", round(x$R,3), "\n",
                    "SE = ", round(x$se,3), "\n",
                    "CI = [", round(x$CI.R[1],3), ", ", round(x$CI.R[2],3), "]", "\n",
                    "P  = ", signif(x$P, 3), "\n\n", 
                    sep="")
                cat("Bootstrapping and Permutation test:", "\n")
                print(format(rbind(x$boot, x$permut), digits = 3))  
        } 
        
        if(x$datatype=="Gaussian" & x$method == "ANOVA") {
                cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n",
                    "Call = ", deparse(x$call), 
                    "\n", "Data: ", x$nobs, " observations and ", x$ngroups, " groups", "\n\n",
                    sep = "")   
                cat("Repeatability: ", "\n", 
                    "R  = ", round(x$R,3), "\n",
                    "SE = ", round(x$se,3), "\n",
                    "CI = [", round(x$CI.R[1],3), ", ", round(x$CI.R[2],3), "]", "\n",
                    "P  = ", signif(x$P[1], 3), " [", attr(x$P, "names")[1], "]", "\n", 
                    "     ", signif(x$P[2], 3), " [", attr(x$P, "names")[2], "]", "\n\n", 
                    sep="")
                cat("Permutation test:", "\n")
                print(format(rbind(x$boot, x$permut), digits = 3))  
        } 
        
        if(x$datatype=="Gaussian" & x$method == "LMM.reml") { 
                cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n",
                    "Call = ", deparse(x$call), 
                    "\n", "Data: ", x$nobs, " observations and ", x$ngroups, " groups", "\n\n",
                    sep = "")   
                cat("Repeatability: ", "\n", 
                    "R  = ", round(x$R,3), "\n",
                    "SE = ", round(x$se,3), "\n",
                    "CI = [", round(x$CI.R[1],3), ", ", round(x$CI.R[2],3), "]", "\n",
                    "P  = ", signif(x$P[1], 3), " [", attr(x$P, "names")[1], "]", "\n", 
                    "     ", signif(x$P[2], 3), " [", attr(x$P, "names")[2], "]", "\n\n", 
                    sep="")
                cat("Bootstrapping and Permutation test:", "\n")
                print(format(rbind(x$boot, x$permut), digits = 3))  
                
                cat("\n", "Likelihood ratio test: ", "\n",
                    "logLik full model = ", x$LRT["LRT.mod"], "\n",
                    "logLik red. model = ", x$LRT["LRT.red"], "\n",
                    "D  = ", signif(x$LRT["LRT.D"], 3), ", " , "df = ", x$LRT["LRT.df"], ", ",
                    "P  = ", signif(x$LRT["LRT.P"], 3), 
                    sep = "")
        } 
	
}	
