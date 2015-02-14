#' Summary of a rpt object
#' 
#' 
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
#' summary(rpt.Weight <- rpt.mcmcLMM(Weight, BirdID))
#' detach(BodySize)
#'       
#' @keywords models
#' 
#' @export
#' 
#' 
#' 
#' 
#' 


summary.rpt <- function(x) {
        #rpt.corr and  rpt.remlLMM and rpt.aov
        if(x$datatype =="Gaussian" & ((x$method == "corr") | (x$method == "LMM.REML"))) {
                # bootstrap and permutation table 
                CI.perm  <- quantile(x$R.permut, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
                x$rpt    <- structure(data.frame(x$R, x$se ,unname(x$P[1]), x$CI.R[1], x$CI.R[2]), 
                                      names = c("R", "SE", attr(x$P, "names")[1], 
                                                attr(CI.perm, "names")[1], attr(CI.perm, "names")[2]),
                                      row.names = "rpt")
                bootperm <- structure(data.frame(c(length(x$R.boot), length(x$R.permut)),
                                               c(mean(x$R.boot), mean(x$R.permut)),
                                               c(median(x$R.boot), median(x$R.permut)),
                                               c(unname(x$CI.R[1]), unname(CI.perm[1])),
                                               c(unname(x$CI.R[2]), unname(CI.perm[2]))),
                                    names = c("N", "Mean", "Median", 
                                              attr(CI.perm, "names")[1], attr(CI.perm, "names")[2]),
                                    row.names = c("boot", "permut"))
                x$boot   <-  bootperm[1, ]
                x$permut <-  bootperm[2, ]
                class(x) <- "summary.rpt"
                return(x) 		
        } 
        
        if(x$datatype=="Gaussian" & x$method == "ANOVA") {
                # anova repeatability and permutation table 
                CI.perm  <- quantile(x$R.permut, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
                x$rpt    <- structure(data.frame(x$R, x$se, unname(x$P[1]), x$CI.R[1], x$CI.R[2]), 
                                      names = c("R", "SE", attr(x$P, "names")[1], 
                                                attr(CI.perm, "names")[1], attr(CI.perm, "names")[2]))
                x$permut <- structure(data.frame(length(x$R.permut), mean(x$R.permut),
                                                 median(x$R.permut), unname(x$P[2]),
                                                 unname(CI.perm[1]), unname(CI.perm[2])),
                                      names = c("N", "Mean", "Median",
                                                attr(x$P, "names")[2], 
                                                attr(CI.perm, "names")[1],
                                                attr(CI.perm, "names")[2]))
                class(x) <- "summary.rpt"
                return(x)         	
        } 
        
        if(x$datatype=="Gaussian" & length(x$P)>1 & length(x$R)>1) {
                warning("Not yet implemented")
#                 cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n")
#                 for(i in 1: length(x$R)) {
#                         cat("Repeatability for ", names(x$R)[i], "\n",
#                             "R  = ", round(x$R[i],3), "\n",
#                             "SE = ", round(x$se[i],3), "\n",
#                             "CI = [", round(x$CI.R[i,1],3), ", ", round(x$CI.R[i,2],3), "]", "\n",
#                             "P  = ", signif(x$P[i,1], 3), " [", attr(x$P, "names")[1], "]", "\n", 
#                             "     ", signif(x$P[i,2], 3), " [", attr(x$P, "names")[2], "]", "\n\n", 
#                             sep="")
#                 }
        }
        
        if(x$datatype!="Gaussian" & x$method=="PQL") {
                # CI for bootstrap and permutation
                CI.boot       <- t(cbind(x$CI.link, x$CI.org))
                CI.perm       <- t(cbind(quantile(x$R.permut$R.link, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE),
                                         quantile(x$R.permut$R.org, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)))         
                # bootstrap table
                Boot.df.temp  <- data.frame(N = c(length(x$R.boot$R.link), length(x$R.boot$R.org)),
                                          mean = c(mean(x$R.boot$R.link), mean(x$R.boot$R.org)), 
                                          median = c(median(x$R.boot$R.link), median(x$R.boot$R.org)),
                                          row.names = c("link", "org"))
                Boot.df       <- cbind(Boot.df.temp, CI.boot)
                # permutation table
                Perm.df.temp  <- data.frame(N = c(length(x$R.permut$R.link), length(x$R.permut$R.org)),
                                            mean = c(mean(x$R.permut$R.link), mean(x$R.permut$R.org)), 
                                            median = c(median(x$R.permut$R.link), median(x$R.permut$R.org)),
                                            row.names = c("link", "org"))
                Perm.df       <- cbind(Perm.df.temp, CI.perm)
                
                x$boot   <- Boot.df
                x$permut <- Perm.df
                class(x) <- "summary.rpt"
                return(x) 
        }
        	
}






