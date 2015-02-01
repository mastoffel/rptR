#' Plot a rpt object
#' 
#' Plots the distribution of repeatability estimates with a confidence interval
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
#' # repeatability estimation for tarsus length - a very high R
#' data(BodySize)
#' attach(BodySize)
#' (rpt.BS <- rpt.remlLMM(Tarsus, BirdID, nboot=100, npermut=100))   
#' # reduced number of nboot and npermut iterations
#' plot(rpt.BS)
#' detach(BodySize)
#'       
#' @keywords models
#' 
#' @export
#' 
#' 
#' 
plot.rpt <- function(x, y, ...) {
        # more default arguments !
        # rpt.corr, rpt.rmlLMM
        if(x$datatype=="Gaussian" & length(x$R)==1) {
                
        CI.l <- unname(x$CI.R[1])
        CI.u <- unname(x$CI.R[2])
        
        boot.plot <- hist(x$R.boot, breaks = "FD",
                     plot = FALSE)
        
        v.pos <- max(boot.plot$counts) # 1/10 * max(boot.plot$counts)
        
        hist(x$R.boot, breaks = "FD", 
             ylim = c(0, v.pos*1.5),
             xlab = "ICC estimates",
             main="Distribution of ICC estimates from bootstrap")
       
        arrows(CI.l, v.pos*1.15, CI.u, v.pos*1.15, length=0.3, angle=90, code=3,
               lwd = 2.5, col = "black")
        lines(x = c(x$R, x$R), y = c(0, v.pos * 1.15),
              lwd = 2.5, col = "grey", lty = 5)
        points(x$R, v.pos*1.15, cex = 1.2, pch = 19, col = "red")
        legend("topleft", pch = 19, cex = 1, bty = "n",
               col = c("red"), c("ICC with CI"),
               box.lty = 0)
        }
        
}
        
#         if(x$datatype=="Gaussian" & length(x$P)==1 & length(x$R)==1) {
#                 cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n",
#                     "R  = ", round(x$R,3), "\n",
#                     "SE = ", round(x$se,3), "\n",
#                     "CI = [", round(x$CI.R[1],3), ", ", round(x$CI.R[2],3), "]", "\n",
#                     "P  = ", signif(x$P, 3), "\n\n", 
#                     sep="")  		
#         } 
#         if(x$datatype=="Gaussian" & length(x$P)>1 & length(x$R)==1) {
#                 cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n",
#                     "R  = ", round(x$R,3), "\n",
#                     "SE = ", round(x$se,3), "\n",
#                     "CI = [", round(x$CI.R[1],3), ", ", round(x$CI.R[2],3), "]", "\n",
#                     "P  = ", signif(x$P[1], 3), " [", attr(x$P, "names")[1], "]", "\n", 
#                     "     ", signif(x$P[2], 3), " [", attr(x$P, "names")[2], "]", "\n\n", 
#                     sep="")
#         } 
#         if(x$datatype=="Gaussian" & length(x$P)==1 & length(x$R)>1) {
#                 print("h3")
#                 cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n")
#                 for(i in 1: length(x$R)) {
#                         cat("Repeatability for ", names(x$R)[i], "\n",
#                             "R  = ", round(x$R[i],3), "\n",
#                             "SE = ", round(x$se[i],3), "\n",
#                             "CI = [", round(x$CI.R[i,1],3), ", ", round(x$CI.R[i,2],3), "]", "\n",
#                             "P  = ", signif(x$P[i], 3), "\n\n", 
#                             sep="")
#                 }
#         }
#         if(x$datatype=="Gaussian" & length(x$P)>1 & length(x$R)>1) {
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
#         }
#         
#         if(x$datatype!="Gaussian" & x$method=="PQL") {
#                 cat("\n", "Repeatability calculation using the ", x$method, " method and ", x$link, "link", "\n\n",
#                     "Esimated overdistpersion (omega) = ", x$omega, "\n\n",
#                     "Link scale repeatabilities:","\n",
#                     "R  = ", round(x$R.link,3), "\n",
#                     "SE = ", round(x$se.link,3), "\n",
#                     "CI = [", round(x$CI.link[1],3), ", ", round(x$CI.link[2],3), "]", "\n",
#                     "P  = ", signif(x$P.link, 3), "\n\n", 
#                     "Original scale repeatabilities:","\n",
#                     "R  = ", round(x$R.org,3), "\n",
#                     "SE = ", round(x$se.org,3), "\n",
#                     "CI = [", round(x$CI.org[1],3), ", ", round(x$CI.org[2],3), "]", "\n",
#                     "P  = ", signif(x$P.org, 3), "\n\n", 
#                     sep="")  
#         }
#         if(x$datatype!="Gaussian" & x$method=="MCMC") {
#                 cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n",
#                     "Link scale repeatabilities:","\n",
#                     "R  = ", round(x$R.link,3), "\n",
#                     "SE = ", round(x$se.link,3), "\n",
#                     "CI = [", round(x$CI.link[1],3), ", ", round(x$CI.link[2],3), "]", "\n",
#                     "P  = ", signif(x$P.link, 3), "\n\n", 
#                     "Original scale repeatabilities:","\n",
#                     "R  = ", round(x$R.org,3), "\n",
#                     "SE = ", round(x$se.org,3), "\n",
#                     "CI = [", round(x$CI.org[1],3), ", ", round(x$CI.org[2],3), "]", "\n",
#                     "P  = ", signif(x$P.org, 3), "\n\n", 
#                     sep="")  
#         }	
	
