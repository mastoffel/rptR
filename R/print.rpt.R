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
        
if (x$ratio == TRUE) {
        header_gaussian <- "Repeatability calculation using the glmm method"
        header_nongaussian <- "Repeatability calculation using the glmm method and"
        header2 <- "Repeatability for "
        PE <- "R  = "
        subheader_link <- "Link scale repeatabilities:"
        subheader_org <- "Original scale repeatabilities:"
} else if (x$ratio == FALSE) {
        header_gaussian <- "Variance estimation using the glmm method"
        header_nongaussian <- "Variance estimation using the glmm method and"
        header2 <- "Estimated Variance for "
        PE <- "Var  = "
        subheader_link <- "Link scale variance:"
        subheader_org <- "Original scale variance:"
}

        
        if (x$datatype == "Gaussian") {
                cat("\n\n")
                cat(header_gaussian, "\n\n")
                for (i in 1:length(x$R)) {
                        cat(header2, names(x$R)[i], "\n",
                                PE, round(unlist(x$R[i]), 3), "\n", sep = "")
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
            cat("\n\n")
            cat(header_nongaussian, x$link, 
                    "link", "\n\n") 
            # grnames <- names(x$R)
            for (i in 1:ncol(x$R)) {
                    grname <- names(x$R)[i]
                    cat(header2, names(x$R)[i], "\n",
                         "--------------------------------", "\n",
                         subheader_link, 
                         "\n", PE, round(x$R["R_link", grname], 3), "\n", "SE = ", round(x$se["se_link", grname], 3),
                          "\n", "CI = [", round(x$CI_emp$CI_link[grname, 1], 3), ", ",  round(x$CI_emp$CI_link[grname, 2], 3), "]", "\n", 
                          "P  = ", signif(x$P[grname,  "P_permut_link"], 3), " [", "Permutation", "]", "\n",
                        "     ", signif(x$P[grname,  "LRT_P"], 3), " [", "LRT", "]", "\n\n",               
                        
                        subheader_org,
                        "\n", PE, round(x$R["R_org", grname], 3), "\n", "SE = ", round(x$se["se_org", grname], 3), "\n", 
                        "CI = [", round(x$CI_emp$CI_org[grname, 1], 3), ", ", round(x$CI_emp$CI_org[grname, 2], 3), "]", "\n", 
                        "P  = ", signif(x$P[grname,  "P_permut_org"], 3), " [", "Permutation", "]", "\n",
                        "     ", signif(x$P[grname,  "LRT_P"], 3), " [", "LRT", "]", "\n\n", sep = "")
                    
    }
    }
}

        
        
        
 
