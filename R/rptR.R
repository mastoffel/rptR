#' rptR: Repeatability calculation for Gaussian and Non-Gaussian data
#'
#' 
#' @description A collection of functions for caluculating point estimates, interval estimates and 
#' significance tests of the repeatability (intra-class correlation coefficient) of measurements.
#' The function \link{rpt} is a the core functions that calls more specialised functions as required. 
#' Specialised functions can also be called directly (see \link{rpt} for details).
#' All functions return lists of values. The function \link{summary.rpt} produces summaries in a detailed format.
#' 
#' @note Currently there is are a number of functions that use the three arguments, (1) \code{data}, specifies
#' the \code{data.frame} containing the variables, \code{y} specifies the response variable and \code{groups} specifies
#' the group variable in the \code{data.frame}. These are used to 
#' estimate standard repeatabilities (see e.g. function \link{rpt}). There is another group of 
#' functions (e.g. \link{rpt.adj}) that uses the arguments \code{formula} and \code{grname} besides the \code{data} argument, which allows 
#' to estimate adjusted repeatabilities (controlling for fixed effects) and the estimation of multiple variance 
#' components simulatneously (multiple random effects). In the long run, the two groups of functions will be merged 
#' and will use the more flexible fomula arguments. So far, adjusted repeatabilities are only implemented for 
#' Gaussian data using REML estimation.
#' 
#' @author Holger Schielzeth (holger.schielzeth@@ebc.uu.se), Shinichi Nakagawa (shinichi.nakagawa@@otago.ac.nz),
#'         Martin Stoffel (martin.adam.stoffel@@gmail.com)
#' 
#' @references Nakagawa, S. and Schielzeth, H. (2010) \emph{Repeatability for Gaussian and non-Gaussian data: 
#'             a practical guide for biologists}. Biological Reviews 85: 935-956
#'             
#' @docType package
#' @name rptR
NULL 
