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
        if(x$datatype =="Gaussian" & ((x$method == "corr") | (x$method == "LMM.reml"))) {
                # bootstrap and permutation table 
                CI.R     <- x$CI.R
                CI.perm  <- quantile(x$R.permut, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
                CI.df    <- t(data.frame(CI.R, CI.perm)) #round(t(data.frame(CI.R, CI.perm)), 4)
                RandBoot.df <- data.frame(N = c(length(x$R.boot), length(x$R.permut)),
                                          mean = c(mean(x$R.boot), mean(x$R.permut)), 
                                          median = c(median(x$R.boot), median(x$R.permut)))
                RandBoot.df <- cbind(RandBoot.df, CI.df)
                row.names(RandBoot.df) <- c("boot", "permut")
                x$boot   <- RandBoot.df[1, ]
                x$permut <- RandBoot.df[2, ]
                class(x) <- "summary.rpt"
                return(x) 		
        } 
        
        if(x$datatype=="Gaussian" & x$method == "ANOVA") {
                # bootstrap and permutation table 
                CI.perm  <- quantile(x$R.permut, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
                x$permut <- cbind(data.frame(N = length(x$R.permut),
                                       mean = mean(x$R.permut), 
                                       median = median(x$R.permut),
                                       row.names = "permut"),
                                       t(CI.perm))
                class(x) <- "summary.rpt"
                return(x)         	
        } 

        	
}






