#' Plot a rpt object
#' 
#' Plots the distribution of repeatability estimates from bootstrapping and permutation tests.
#' 
#' @param x An rpt object returned from one of the rpt functions.
#' @param grname The name of the grouping factor to plot.
#' @param scale Either "link" or "original" scale results for results of non-Gaussian functions.
#' @param type Either "boot" or "permut" for plotting the results of bootstraps or permutations.
#' @param main Plot title
#' @param breaks hist() argument
#' @param xlab x-axis title
#' @param \dots Additional arguments to the hist() function for customized plotting.
#'
#' @references 
#' Nakagawa, S. & Schielzeth, H. (2010) \emph{Repeatability for Gaussian and 
#' non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
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
plot.rpt <- function(x, grname = names(x$ngroups), scale = c("link", "original"), type = c("boot", "permut"), 
        main = NULL, breaks = "FD", xlab = NULL, ...) {
    
    # save ellipsis args
    dots <- list(...)
    
    if (is.null(xlab)) {
            if (x$ratio==TRUE){  
                    xlab = "Repeatability estimates" 
            } else {
                    xlab = "Variance estimates"
            }
    }
        
    # initialising
    if (length(type) != 1) type <- type[1]
    if (length(scale) != 1) scale <- scale[1]
    if (length(grname) != 1) grname <- grname[1]
    
    if (x$datatype != "Gaussian") {
        if (is.null(main)) {
            if (type == "boot") {
                if (scale == "link") {
                        if (x$ratio == FALSE){
                                main <- paste("Link-scale approximation bootstrap \nvariance estimates for", grname)
                        } else {
                                main <- paste("Link-scale approximation bootstrap \nrepeatabilities for", grname)
                        }
                }
                if (scale == "original") {
                        if (x$ratio == FALSE) stop("No original-scale approximation for Variances. Change to scale = 'link'")
                         main <- paste("Original-scale approximation \nrepeatabilities for", grname)
                }
            } else if (type == "permut") {
                    if (grname == "Residual" | grname == "Overdispersion") stop("No permutation tests for Residual or Overdispersion")
                if (scale == "link") {
                        main <- paste("Link-scale approximation permutation test \nrepeatabilities for", grname)   
                }
                if (scale == "original") {
                        if (x$ratio == FALSE) stop("No original-scale approximation for Variances. Change to scale = 'link'")
                        main <- paste("Original-scale approximation permutation test \nrepeatabilities for", grname)
                }
                 
                
            }
        }
    }
    
    
    
    # make bootstrap histogram
    boot_hist <- function(R, R.boot, CI.l, CI.u, xlab. = xlab, breaks. = breaks, main. = main, 
        ...) {
        dots <- list(...)
        
        # y position of confidence band
        v.pos <- max((graphics::hist(R.boot, breaks = breaks., plot = FALSE))$counts)
        # plot
        do.call(graphics::hist, args = c(list(R.boot, breaks = breaks., ylim = c(0, v.pos * 1.5), xlab = xlab., main = main.), dots))
        graphics::lines(x = c(R, R), y = c(0, v.pos * 1.15), lwd = 1.5, col = "grey", lty = 5)
        graphics::arrows(CI.l, v.pos * 1.15, CI.u, v.pos * 1.15, length = 0.05, angle = 90, code = 3, 
            lwd = 1.5, col = "black")
        graphics::points(R, v.pos * 1.15, cex = 1.1, pch = 19, col = "cornflowerblue")
        if(x$ratio == FALSE){
                graphics::legend("topleft", pch = 19, cex = 0.8, bty = "n", col = c("cornflowerblue"), c("Variance estimate with CI"), 
                        box.lty = 0)      
        } else {
                graphics::legend("topleft", pch = 19, cex = 0.8, bty = "n", col = c("cornflowerblue"), c("Repeatability with CI"), 
                        box.lty = 0)
        }
    
    }
    
    permut_hist <- function(R, R.permut, xlab. = xlab, CI = x$CI, breaks. = breaks, main. = main, 
        ...) {
        dots <- list(...)
        # get CI for permutation
        CI.perm <- stats::quantile(R.permut, c((1 - CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
        Median.R <- stats::median(R.permut)
        # y position of confidence band
        v.pos <- max((graphics::hist(R.permut, breaks = breaks., plot = FALSE))$counts)
        # plot
        do.call(graphics::hist, args = c(list(R.permut, breaks = breaks., ylim = c(0, v.pos * 1.5), xlab = xlab., main = main.), dots))
        graphics::lines(x = c(Median.R, Median.R), y = c(0, v.pos * 1.15), lwd = 1.5, col = "grey", 
            lty = 5)
        graphics::lines(x = c(R, R), y = c(0, v.pos * 1.3), lwd = 1.5, col = "grey", lty = 5)
        graphics::arrows(unname(CI.perm[1]), v.pos * 1.15, unname(CI.perm[2]), v.pos * 1.15, length = 0.05, 
            angle = 90, code = 3, lwd = 1.5, col = "black")
        graphics::points(Median.R, v.pos * 1.15, cex = 1.2, pch = 19, col = "black")
        graphics::points(R, v.pos * 1.3, cex = 1.1, pch = 19, col = "cornflowerblue")
        graphics::legend("topleft", pch = 19, cex = 0.8, bty = "n", col = c("black", "cornflowerblue"), c("Median of repeatabilities from permuted datasets with 95% percentiles", 
                "Observed repeatability"), box.lty = 0)
    }
    
    
    if (x$datatype == "Poisson" | x$datatype == "Binary" | x$datatype == "Proportion") {
            if (type == "boot") {
                    if (scale == "link") {
                            boot_hist(R = x$R[2, grname], R.boot = unname(unlist(x$R_boot_link[grname])), 
                                    CI.l = unname(x$CI_emp$CI_link[grname, 1]), 
                                    CI.u = unname(x$CI_emp$CI_link[grname, 2]), 
                                    main. = main, ...)   
                    } else if (scale == "original") {
                            boot_hist(R = x$R[2, grname], R.boot = unname(unlist(x$R_boot_org[grname])), 
                                    CI.l = unname(x$CI_emp$CI_org[grname, 1]), 
                                    CI.u = unname(x$CI_emp$CI_org[grname, 2]), 
                                    main. = main, ...)   
                    }
            } else if (type == "permut") {
                    if (scale == "link") {
                            permut_hist(R = x$R[2, grname], R.permut = unname(unlist(x$R_permut_link[grname])), 
                                    main. = main, ...)  
                    }
                    if (scale == "org") {
                            permut_hist(R = x$R[2, grname], R.permut = unname(unlist(x$R_permut_org[grname])), 
                                    main. = main, ...)  
                    }
            }
    }
    
    if (x$datatype == "Gaussian") {
            if (is.null(main)) {
                    if (type == "boot") {
                            if (scale == "link") {
                                    if (x$ratio == FALSE){
                                            main <- paste("Bootstrap variance estimates \nfor", grname)   
                                    } else {
                                            main <- paste("Bootstrap repeatabilities \nfor", grname)
                                    }
                            }
             
                    } else if (type == "permut") {
                            if (grname == "Residual" | grname == "Overdispersion") stop("No permutation tests for Residual or Overdispersion")
                            if (x$ratio == FALSE){
                                    main <- paste("Permutation variance estimates \nfor", grname)   
                            } else {
                                    main <- paste("Permutation repeatabilities \nfor", grname)
                            }
                    }
            }      
            
            
       if (type == "boot") {
                boot_hist(R = x$R[grname], R.boot = unlist(x$R_boot[grname]), CI.l = unname(x$CI_emp[grname, 
                  1]), CI.u = unname(x$CI_emp[grname, 2]), main. = main, ...)
            }
       if (type == "permut") {
                permut_hist(R = x$R[grname], R.permut = unlist(x$R_permut[grname]), main. = main, ...)
       }
    }
    
}
 
