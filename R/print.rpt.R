#' Print a rpt object
#' 
#' Displays the results a rpt object (i.e. the result of a rpt function call) in a nice form.
#' 
#' @param x An rpt object returned from one of the rpt functions
#' @param \dots Additional arguments; none are used in this method.
#'
#' @references 
#' Nakagawa, S. and Schielzeth, H. (2010) \emph{Repeatability for Gaussian and 
#'              non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
#' 
#' @author Holger Schielzeth  (holger.schielzeth@@ebc.uu.se), 
#'         Shinichi Nakagawa (shinichi.nakagawa@@otago.ac.nz) &
#'         Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'      
#' @seealso \link{rpt}, \link{rpt.corr}, \link{rpt.aov}, \link{rpt.remlLMM}, \link{rpt.mcmcLMM},
#'          \link{rpt.binomGLMM.add}, \link{rpt.binomGLMM.multi}, \link{rpt.poisGLMM.add}, \link{rpt.poisGLMM.multi}
#' 
#'       
#' @keywords models
#' 
#' @export
#' 
#' 
#' @examples  
#' # repeatability estimation for weight 
#' data(BodySize)
#' rpt.Weight <- rpt.mcmcLMM(BodySize, Weight, BirdID)
#' print(rpt.Weight)  # alternative call to print.rpt() through print()
#' 
print.rpt <- function(x, ...) {
    # rpt.corr, rpt.mcmcLMM

    # rpt.remlLMM.adj   
    if (x$datatype == "Gaussian" & length(x$P) > 1 & length(x$R) > 1) {
        cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n")
        for (i in 1:length(x$R)) {
            cat("Repeatability for ", names(x$R)[i], "\n", "R  = ", round(x$R[i], 3), 
                "\n", "SE = ", round(x$se[i], 3), "\n", "CI = [", round(x$CI.R[i, 1], 
                  3), ", ", round(x$CI.R[i, 2], 3), "]", "\n", "P  = ", signif(x$P[i, 
                  1], 3), " [", attr(x$P, "dimnames")[[2]][1], "]", "\n", "     ", signif(x$P[i, 
                  2], 3), " [", attr(x$P, "dimnames")[[2]][2], "]", "\n\n", sep = "")
        }
    }
        
    if (x$datatype == "Poisson" | x$datatype == "Binary") {
            cat("\n", "Repeatability calculation using the glmm method and", x$link, 
                    "link",  "\n\n",  "Estimated overdispersion = ", x$overdisp, "\n\n") 
            # grnames <- names(x$R)
            for (i in 1:ncol(x$R)) {
                    grname <- names(x$R)[i]
                    cat("Repeatability for ", grname, "\n",
                         "--------------------------------", "\n",
                         "Link scale repeatabilities:", 
                         "\n", "R  = ", round(x$R["R_link", grname], 3), "\n", "SE = ", round(x$se[grname, "se_link"], 3),
                          "\n", "CI = [", round(x$CI_emp$CI_link[grname, 1], 3), ", ",  round(x$CI_emp$CI_link[grname, 2], 3), "]", "\n", 
                          "P  = ", signif(x$P[grname,  "P_permut_link"], 3), "\n\n",
                        "Original scale repeatabilities:",
                        "\n", "R  = ", round(x$R["R_org", grname], 3), "\n", "SE = ", round(x$se[grname, "se_org"], 3), "\n", 
                        "CI = [", round(x$CI_emp$CI_org[grname, 1], 3), ", ", round(x$CI_emp$CI_org[grname, 2], 3), "]", "\n", 
                        "P  = ", signif(x$P[grname,  "P_permut_org"], 3), "\n\n", sep = "")
            }
    }
        

}
 
