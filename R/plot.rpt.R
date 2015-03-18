#' Plot a rpt object
#' 
#' Plots the distribution of repeatability estimates from bootstrapping and permutation tests.
#' For MCMC methods
#' 
#' @param x An rpt object returned from one of the rpt functions.
#' @param type Showing either the bootstap or permutation test results.
#' @param scale Either link or original scale results for binomial or poisson data and the multiplicative overdispersion model.
#' @param main Plot title
#' @param breaks hist() argument
#' @param xlab x-axis title
#' @param \dots Additional arguments to the hist() function.
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
#' (rpt.BS <- rpt.remlLMM("Tarsus", "BirdID", data = BodySize, nboot=100, npermut=100))   
#' # reduced number of nboot and npermut iterations
#' plot(rpt.BS)
#'       
#' @keywords models
#' 
#' @export
#' 
#' 
#' 
plot.rpt <- function(x, type = c("boot", "permut"), scale = c("link", "original"), 
                     main = NULL, breaks = "FD", xlab = "Repeatability estimates", ...) {
        
        # initialising
        if (length(type) != 1)  type  <- type[1]
        if (length(scale) != 1) scale <- scale[1]
        
        if (x$datatype!="Gaussian" & x$method=="PQL") {
                if (is.null(main)){
                        if (type == "boot") {
                                if (scale == "link") main <- "Link scale distribution of repeatability estimates from bootstrap"
                                if (scale == "original") main <- "Original scale distribution of repeatability estimates from bootstrap"
                        } else if (type == "permut") {
                                if (scale == "link") main <- "Link scale distribution of repeatability estimates from permutation"
                                if (scale == "original") main <- "Original scale distribution of repeatability estimates from permutation"
                        
                        }
                }
        }
        if(x$datatype=="Gaussian" & ((x$method == "corr") | (x$method == "LMM.REML") | (x$method == "ANOVA"))) {
                if (is.null(main)){
                        if (type == "boot") main <- "Distribution of repeatability estimates from bootstrap"
                        if (type == "permut") main <- "Distribution of repeatability estimates from permutation"
                }
        }
        
        # make bootstrap histogram
        boot_hist <- function(R, R.boot, CI.l, CI.u, xlab. = xlab, breaks. = breaks, main. = main, ...) {
                # y position of confidence band
                v.pos <- max((hist(R.boot, breaks = breaks., plot = FALSE))$counts)  
                # plot
                hist(R.boot, breaks = breaks., ylim = c(0, v.pos*1.5), xlab = xlab.,
                     main=main.)
                lines(x = c(R, R), y = c(0, v.pos * 1.15), lwd = 2.5, col = "grey", lty = 5)
                arrows(CI.l, v.pos*1.15, CI.u, v.pos*1.15, 
                       length=0.3, angle=90, code=3, lwd = 2.5, col = "black")
                points(R, v.pos*1.15, cex = 1.2, pch = 19, col = "red")
                legend("topleft", pch = 19, cex = 1, bty = "n", col = c("red"), 
                       c("Repeatability with CI"), box.lty = 0)
        }
        
        permut_hist <- function(R, R.permut, xlab. = xlab, CI = x$CI, breaks. = breaks, main. = main, ...) {
                # get CI for permutation
                CI.perm  <- quantile(R.permut, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
                Median.R <- median(R.permut)
                # y position of confidence band
                v.pos <- max((hist(R.permut, breaks = breaks., plot = FALSE))$counts)  
                # plot
                hist(R.permut, breaks = breaks., ylim = c(0, v.pos*1.5), xlab = xlab.,
                     main=main.)
                lines(x = c(Median.R, Median.R), y = c(0, v.pos * 1.15), lwd = 2.5, col = "grey", lty = 5)
                lines(x = c(R, R), y = c(0, v.pos * 1.3), lwd = 2.5, col = "grey", lty = 5)
                arrows(unname(CI.perm[1]), v.pos*1.15, unname(CI.perm[2]), v.pos*1.15, 
                       length=0.3, angle=90, code=3, lwd = 2.5, col = "black")
                points(Median.R, v.pos*1.15, cex = 1.2, pch = 19, col = "black")
                points(R, v.pos*1.3, cex = 1.2, pch = 19, col = "red")
                legend("topleft", pch = 19, cex = 1, bty = "n", col = c("black", "red"), 
                       c("Median of simulated repeatabilities with CI", "Observed repeatability"), box.lty = 0)
        }
        
        if((x$datatype=="Gaussian") & (x$method == "LMM.REML") & (length(x$R)>1)) {
                devAskNewPage(ask = TRUE)
                on.exit(devAskNewPage(ask = NULL))
                if(type == "boot") {
                        for (i in 1:length(x$R)) {
                                boot_hist(R = x$R[i], R.boot = x$R.boot[i, ], 
                                        CI.l = unname(x$CI.R[i, 1]),
                                        CI.u = unname(x$CI.R[i, 2]), 
                                        main. = paste("Bootstrap repeatabilities for", names(x$R)[i]), ...)
                        }
                } else if(type == "permut") {
                        for (i in 1:length(x$R)) {
                                permut_hist(R = x$R[i], R.permut = x$R.permut[i, ], 
                                            main. = paste("Permutation repeatabilities for", names(x$R)[i]), ...)
                        }
                } 
        }
        
        if(x$datatype=="Gaussian" & ((x$method == "corr") | (x$method == "LMM.REML")) & length(x$R)==1) {
                if(type == "boot") {
                        boot_hist(R = x$R, R.boot = x$R.boot, 
                                  CI.l = unname(x$CI.R[1]),
                                  CI.u = unname(x$CI.R[2]), ...)
                } else if(type == "permut") {
                        if (x$method == "corr") {
                                permut_hist(R = x$R, R.permut = x$R.permut, ...) # x$P??
                        } else if (x$method == "LMM.REML") {
                                # no red point. unclear which p to plot
                                permut_hist(R = x$R, R.permut = x$R.permut, ...)      
                        }
                } 
        }
        
        
        if(x$datatype=="Gaussian" & x$method == "ANOVA") {
              if (type == "boot") {
                      warning("type = 'boot' not available. Use type = 'permut'.")
              }
              if (type == "permut") {
              permut_hist(R = x$R, R.permut = x$R.permut, ...)  
              }
        } 
        
        if(x$datatype!="Gaussian" & x$method=="PQL") {
                # if (is.null(scale)) warning("Set scale to 'link' or 'original' to show the respective plot")
                if (scale == "link") {
                        if(type == "boot") {
                                boot_hist(R = x$R.link, R.boot = x$R.boot$R.link, 
                                          CI.l = unname(x$CI.link[1]),
                                          CI.u = unname(x$CI.link[2]), ...)
                        } else if (type == "permut") {
                                        permut_hist(R = x$R.link, R.permut = x$R.permut$R.link, ...)  
                                }
                } 
                if (scale == "original") {
                        if(type == "boot") {
                                boot_hist(R = x$R.org, R.boot = x$R.boot$R.org, 
                                          CI.l = unname(x$CI.link[1]),
                                          CI.u = unname(x$CI.link[2]), ...)
                        } else if (type == "permut") {
                                permut_hist(R = x$R.org, R.permut = x$R.permut$R.org, ...)  
                        }
                }
                
        }
        
        if(x$method=="LMM.MCMC" | x$method=="MCMC") { 
                plot(x$mod)
        }
        
}
        

