#' Summary of a rpt object
#'
#' 
#' @param object An rpt object returned from one of the rpt functions
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
#' @keywords models
#' 
#' @export
#' 
#' 
#' 
#' 
#' 

summary.rpt <- function(object, ...) {

        # helper functions for 
        CI <- R_est$CI
        calc_CI <- function(x) {
                out <- quantile(x, c((1 - CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
                out
        }
        
        extr_comps <- function(x) {
                CI <- calc_CI(x)
                out <- data.frame("N" = length(x), "Mean" = mean(x), "Median" = median(x),
                        calc_CI(x)[[1]], calc_CI(x)[[2]])
                names(out)[4:5] <- names(calc_CI(x))
                out
        }
        
        if(object$datatype=="Gaussian") {
                # warning("Not yet implemented")
#                 cat("\n", "Repeatability calculation using the ", object$method, " method", "\n\n")
        for(i in 1:length(object$R)) {
                CI.perm  <- quantile(object$R.permut[i, ], c((1-object$CI)/2,1-(1-object$CI)/2), na.rm=TRUE) # should be taken out of the loop to the top when impl.
                object$rpt[[i]]    <- structure(data.frame(object$R[i], object$se[i] ,unname(object$P[i, 1]), object$CI.R[i, 1], object$CI.R[i, 2]), 
                                           names = c("R", "SE", colnames(object$P)[1], 
                                                   colnames(object$CI.R)[1], colnames(object$CI.R)[2]),
                                                   row.names = "rpt")
                bootperm      <- structure(data.frame(c(ncol(object$R.boot), ncol(object$R.permut)),
                                                 c(mean(object$R.boot[i, ]), mean(object$R.permut[i, ])),
                                                 c(median(object$R.boot[i, ]), median(object$R.permut[i, ])),
                                                 c(unname(object$CI.R[i, 1]), unname(CI.perm[1])), # should be CI.perm[i, 1] when implemented
                                                 c(unname(object$CI.R[i, 2]), unname(CI.perm[2]))),
                                                 names = c("N", "Mean", "Median", 
                                                 colnames(object$CI.R)[1], colnames(object$CI.R)[2]),
                                                 row.names = c("boot", "permut"))
                object$boot[[i]]   <-  bootperm[1, ]
                object$permut[[i]] <-  bootperm[2, ]
        }
                class(object) <- "summary.rpt"
                return(object)
        }
        
        
        
        if(object$datatype=="Poisson" | object$datatype=="Binary" |  object$datatype=="Proportion") {
                # warning("Not yet implemented")
                #                 cat("\n", "Repeatability calculation using the ", object$method, " method", "\n\n")
                boot <- list(object$R_boot_org, object$R_boot_link)
                perm <- list(object$R_permut_org, object$R_permut_link)
                for(i in 1:ncol(object$R)) {
                        object$rpt[[i]] <- structure(data.frame(R = object$R[[i]], object$se[[i]], 
                                          do.call(rbind, lapply(object$CI_emp, function(x) x[i, ])),
                                          object$P[i, c("P_permut_org", "P_permut_link")]),
                                           names = c("R", "SE", names(calc_CI(object$R_boot_org[[1]])), "P_val"),
                                          row.names = c("Org", "Link"))
                        object$boot[[i]] <- structure(do.call(rbind, 
                                            lapply(boot, function(x) extr_comps(x[[i]]))),
                                            row.names = c("Org", "Link"))
                        object$permut[[i]] <- structure(cbind(
                                do.call(rbind, lapply(perm, function(x) extr_comps(x[[i]]))),
                                object$P[i, c("P_permut_org", "P_permut_link")]),
                                              row.names = c("Org", "Link"))
                        names(object$permut[[i]])[6] <- "P_val"
                       
                }
                class(object) <- "summary.rpt"
                return(object)
        }
        	
}






