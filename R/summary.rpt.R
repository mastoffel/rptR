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
#' 
#' 

summary.rpt <- function(object, ...) {

        # helper functions for 
        CI <- object$CI
        calc_CI <- function(x) {
                out <- stats::quantile(x, c((1 - CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
                out
        }
        
        extr_comps <- function(x) {
                CI <- calc_CI(x)
                out <- data.frame("N" = length(x), "Mean" = mean(x), "Median" = stats::median(x),
                        calc_CI(x)[[1]], calc_CI(x)[[2]])
                names(out)[4:5] <- names(calc_CI(x))
                out
        }
        
        if(object$datatype=="Gaussian") {
                for(i in 1:length(object$R)) {
                      
                        object$rpt[[i]]    <- structure(data.frame(object$R[i], unlist(object$se)[[i]], as.numeric(matrix(object$CI_emp, ncol = 2)[i, 1]), 
                                                        as.numeric(matrix(object$CI_emp, ncol = 2)[i, 2]), 
                                                        unname(object$P[i, 2]), round(unname(object$P[i, 1]), 3)), 
                                                        names = c("R", "SE", names(object$CI_emp)[1], names(object$CI_emp)[2],
                                                        names(object$P)[2],  names(object$P)[1]), 
                                                           row.names = "rpt")
                        bootperm      <- structure(data.frame(do.call(rbind, lapply(list(object$R_boot[[i]], object$R_permut[[i]]), extr_comps))),
                                         row.names = c("boot", "permut"), names = c("N", "Mean", "Median", names(object$CI_emp)))
        
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
                                          t(object$P[i, c("P_permut_org", "P_permut_link")])),
                                           names = c("R", "SE", names(calc_CI(object$R_boot_org[[1]])), "P_permut"),
                                          row.names = c("Org", "Link"))
                        object$boot[[i]] <- structure(do.call(rbind, 
                                            lapply(boot, function(x) extr_comps(x[[i]]))),
                                            row.names = c("Org", "Link"))
                        object$permut[[i]] <- structure(cbind(
                                do.call(rbind, lapply(perm, function(x) extr_comps(x[[i]]))),
                                t(object$P[i, c("P_permut_org", "P_permut_link")])),
                                              row.names = c("Org", "Link"))
                        names(object$permut[[i]])[6] <- "P_permut"
                       
                }
                class(object) <- "summary.rpt"
                return(object)
        }
        	
}






