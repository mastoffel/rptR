#' Summary of a rpt object
#' 
#' 
#' 
#' @param object An rpt object returned from one of the rpt functions
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

summary.rpt <- function(object, ...) {
        #rpt.corr and  rpt.remlLMM and rpt.aov, rpt.remlLMM.adj with one groups
        if(object$datatype =="Gaussian" & ((object$method == "corr") | (object$method == "LMM.REML")) & length(object$R)==1) {
                # bootstrap and permutation table 
                CI.perm  <- quantile(object$R.permut, c((1-object$CI)/2,1-(1-object$CI)/2), na.rm=TRUE)
                object$rpt    <- structure(data.frame(object$R, object$se ,unname(object$P[1]), object$CI.R[1], object$CI.R[2]), 
                                      names = c("R", "SE", attr(object$P, "names")[1], 
                                                attr(CI.perm, "names")[1], attr(CI.perm, "names")[2]),
                                      row.names = "rpt")
                bootperm <- structure(data.frame(c(length(object$R.boot), length(object$R.permut)),
                                               c(mean(object$R.boot), mean(object$R.permut)),
                                               c(median(object$R.boot), median(object$R.permut)),
                                               c(unname(object$CI.R[1]), unname(CI.perm[1])),
                                               c(unname(object$CI.R[2]), unname(CI.perm[2]))),
                                    names = c("N", "Mean", "Median", 
                                              attr(CI.perm, "names")[1], attr(CI.perm, "names")[2]),
                                    row.names = c("boot", "permut"))
                object$boot   <-  bootperm[1, ]
                object$permut <-  bootperm[2, ]
                class(object) <- "summary.rpt"
                return(object) 		
        } 
        
        if(object$datatype=="Gaussian" & object$method == "ANOVA") {
                # anova repeatability and permutation table 
                CI.perm  <- quantile(object$R.permut, c((1-object$CI)/2,1-(1-object$CI)/2), na.rm=TRUE)
                object$rpt    <- structure(data.frame(object$R, object$se, unname(object$P[1]), object$CI.R[1], object$CI.R[2]), 
                                      names = c("R", "SE", attr(object$P, "names")[1], 
                                                attr(CI.perm, "names")[1], attr(CI.perm, "names")[2]))
                object$permut <- structure(data.frame(length(object$R.permut), mean(object$R.permut),
                                                 median(object$R.permut), unname(object$P[2]),
                                                 unname(CI.perm[1]), unname(CI.perm[2])),
                                      names = c("N", "Mean", "Median",
                                                attr(object$P, "names")[2], 
                                                attr(CI.perm, "names")[1],
                                                attr(CI.perm, "names")[2]))
                object$anovatab <- structure(as.data.frame(object$mod), row.names = c(as.character(object$call)[3], "Residuals"))
                class(object) <- "summary.rpt"
                return(object)         	
        } 
        
        if(object$datatype=="Gaussian" & object$method == "LMM.MCMC") { 
                # Rpt table
                object$rpt    <- structure(data.frame(object$R, object$se, unname(object$P[1]), object$CI.R[1], object$CI.R[2]), 
                                      names = c("R", "SE", "P", 
                                                attr(object$CI.R, "names")[1], attr(object$CI.R, "names")[2]))
                class(object) <- "summary.rpt"
                return(object)  
        } 
        
        if(object$datatype=="Gaussian" & length(object$P)>1 & length(object$R)>1) {
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
        
     
        if(object$datatype!="Gaussian" & object$method=="PQL") {
                # CI for permutation and bootstrap
                CI.perm       <- as.data.frame(rbind(quantile(object$R.permut$R.link, c((1-object$CI)/2,1-(1-object$CI)/2), na.rm=TRUE),
                                         quantile(object$R.permut$R.org, c((1-object$CI)/2,1-(1-object$CI)/2), na.rm=TRUE)))
                CI.boot       <- as.data.frame(rbind(object$CI.link, object$CI.org))
                # Rpt table with link and original scale
                object$rpt    <- structure(data.frame(c(object$R.link,  object$R.org), 
                                                 c(object$se.link, object$se.org),
                                                 c(object$P.link, object$P.org),
                                                 c(object$CI.link[1], object$CI.org[1]),
                                                 c(object$CI.link[2], object$CI.org[2])),
                                      names = c("R", "SE", "P", 
                                                attr(CI.perm, "names")[1], attr(CI.perm, "names")[2]),
                                      row.names = c("link", "original"))
                # CI for bootstrap and permutation
                   
                # bootstrap table
                object$boot <- structure(data.frame(c(length(object$R.boot$R.link), length(object$R.boot$R.org)),
                                             c(mean(object$R.boot$R.link), mean(object$R.boot$R.org)), 
                                             c(median(object$R.boot$R.link), median(object$R.boot$R.org)),
                                             CI.boot[1], CI.boot[2]),
                                             row.names = c("link", "original"),
                                             names = c("N", "Mean", "Median", names(CI.boot)[1], names(CI.boot)[2]))
                # permutation table
                object$permut  <- structure(data.frame(c(length(object$R.permut$R.link), length(object$R.permut$R.org)),
                                                      c(mean(object$R.permut$R.link), mean(object$R.permut$R.org)), 
                                                      c(median(object$R.permut$R.link), median(object$R.permut$R.org)),
                                                      CI.perm[1], CI.perm[2]),
                                           row.names = c("link", "original"),
                                           names = c("N", "Mean", "Median", names(CI.perm)[1], names(CI.perm)[2]))
                class(object) <- "summary.rpt"
                return(object) 
        }
        
        if(object$datatype!="Gaussian" & object$method=="MCMC") { #object$datatype!="Gaussian" & 
                object$rpt    <- structure(data.frame(c(object$R.link,  object$R.org), 
                                                 c(object$se.link, object$se.org),
                                                 c(object$P.link, object$P.org),
                                                 c(object$CI.link[1], object$CI.org[1]),
                                                 c(object$CI.link[2], object$CI.org[2])),
                                      names = c("R", "SE", "P", 
                                                attr(CI.link, "names")[1], attr(CI.link, "names")[2]),
                                      row.names = c("link", "original"))
                class(object) <- "summary.rpt"
                return(object) 
        }
        	
}






