#' Print a rpt object
#' 
#' Displays the results a rpt object (i.e. the result of a rpt function call) in a nice form.
#' 
#' @param x An rpt object returned from one of the rpt functions
#' @param \dots Additional arguments; none are used in this method.
#'
#' @references 
#' Nakagawa, S. & Schielzeth, H. (2010) \emph{Repeatability for Gaussian and 
#'              non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
#' 
#' @author Holger Schielzeth  (holger.schielzeth@@uni-jena.de), 
#'         Shinichi Nakagawa (s.nakagawa@unsw.edu.au),
#'         Martin Stoffel (martin.adam.stoffel@@gmail.com) 
#'      
#'       
#' @keywords models
#' 
#' @export
#' 
#' 
#' 
print.rpt <- function(x, ...) {
    # rpt.corr, rpt.mcmcLMM
#         if (x$datatype == "Gaussian") {
#                         #  rpt.remlLMM.adj with one group factor
#                         cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n", 
#                                 "R  = ", round(x$R, 3), "\n", "SE = ", round(x$se, 3), "\n", "CI = [", round(x$CI.R[1], 
#                                         3), ", ", round(x$CI.R[2], 3), "]", "\n", "P  = ", signif(x$P[1], 3), 
#                                 " [",names(x$P[1]), "]", "\n", "     ", signif(x$P[2], 3), " [", 
#                                 x$P[2], "]", "\n\n", sep = "")     
#                 }
    
    # if (x$datatype == "Gaussian") {
    #     cat("\n", "Repeatability calculation using the glmm method", "\n\n")
    #     for (i in 1:length(x$R)) {
    #         cat("Repeatability for ", names(x$R)[i], "\n",
    #              "R  = ", round(unlist(x$R[i]), 3), "\n",
    #              "SE = ", round(unlist(x$se[i, ]), 3), "\n",
    #              "CI = [", round(x$CI_emp[i, 1], 3), ", ", round(x$CI_emp[i, 2], 3), "]", "\n",
    #               "P  = ", signif(x$P[i, 2], 3), " [", "Permutation", "]", "\n",
    #               "     ", signif(x$P[i, 1], 3), " [", "LRT", "]", "\n\n", sep = "")
    #     }
    # }
        
        
        if (x$datatype == "Gaussian") {
                cat("\n", "Repeatability calculation using the glmm method", "\n\n")
                for (i in 1:length(x$R)) {
                        cat("Repeatability for ", names(x$R)[i], "\n",
                                "R  = ", round(unlist(x$R[i]), 3), "\n", sep = "")
                        if (is.na(unlist(x$se)[1])){ # check if permutations have been done
                        cat("SE = ", x$se, "\n")
                        } else {
                        cat("SE = ", round(unlist(x$se[i, ]), 3), "\n", sep = "")    
                        }
                        if (is.na(unlist(x$CI_emp)[1])){ # check if bootstraps have been done
                        cat("CI = [", x$CI_emp[1], ", ", x$CI_emp[2], "]", "\n", sep = "")
                        } else {
                        cat("CI = [", round(x$CI_emp[i, 1], 3), ", ", round(x$CI_emp[i, 2], 3), "]", "\n", sep = "")
                        }
                        cat("P  = ", signif(x$P[i, 2], 3), " [", "Permutation", "]", "\n",
                        "     ", signif(x$P[i, 1], 3), " [", "LRT", "]", "\n\n", sep = "")
                }
        }
        
        
        
    if (x$datatype == "Poisson" | x$datatype == "Binary" | x$datatype == "Proportion") {
            cat("\n", "Repeatability calculation using the glmm method and", x$link, 
                    "link", "\n\n") 
            # grnames <- names(x$R)
            for (i in 1:ncol(x$R)) {
                    grname <- names(x$R)[i]
                    cat("Repeatability for ", names(x$R)[i], "\n",
                         "--------------------------------", "\n",
                         "Link scale repeatabilities:", 
                         "\n", "R  = ", round(x$R["R_link", grname], 3), "\n", "SE = ", round(x$se["se_link", grname], 3),
                          "\n", "CI = [", round(x$CI_emp$CI_link[grname, 1], 3), ", ",  round(x$CI_emp$CI_link[grname, 2], 3), "]", "\n", 
                          "P  = ", signif(x$P[grname,  "P_permut_link"], 3), " [", "Permutation", "]", "\n\n", 
                    
                        "Original scale repeatabilities:",
                        "\n", "R  = ", round(x$R["R_org", grname], 3), "\n", "SE = ", round(x$se["se_org", grname], 3), "\n", 
                        "CI = [", round(x$CI_emp$CI_org[grname, 1], 3), ", ", round(x$CI_emp$CI_org[grname, 2], 3), "]", "\n", 
                        "P  = ", signif(x$P[grname,  "P_permut_org"], 3), " [", "Permutation", "]", "\n\n",
                    
                        "Likelihood ratio test: P = ", signif(x$P[grname,  "LRT_P"], 3),  "\n\n", sep = "")
            }
    }
        
        

}
 
