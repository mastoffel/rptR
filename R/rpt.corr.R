#' Correlation-based Repeatability
#' 
#' Calculates repeatability based on inter-class correlations
#' 
#' @param y Vector of measurements. Missing values are not allowed.
#' @param groups Vector of group identitities (will be converted to a factor).
#'        Note that each group identity has to appear exactly twice.
#' @param CI Width of the confidence interval between 0 and 1 (defaults to 0.95).
#' @param nboot Number of bootstrapping runs used when calculating an asymptotic
#'        confidence interval (defaults to 1000).
#' @param npermut Number of permutations used when calculating asymptotic 
#'        \emph{P} values (defaults to 1000).
#'   
#' 
#' @return Returns an object of class rpt that is a a list with the following elements: 
#' \item{datatype}{Response distribution (here: "Gaussian").}
#' \item{method}{Method used to calculate repeatability (here: "corr").}
#' \item{R}{Point estimate for repeatability \emph{R}.}
#' \item{se}{Asymptotic standard error for repeatability based on non-parametric bootstrapping.}
#' \item{CI}{Asymptotic confidence interval for repeatability based on non-parametric bootstrapping.}
#' \item{P}{Asymptotic \emph{P} value from a significance test for the intraclass correlation based on permutation.}
#' \item{R.permut}{Permutation samples for \emph{R}.}
#'
#'
#' @references 
#' Sokal, R. R. and F. J. Rohlf (1995). \emph{Biometry: The principles and practice of statistics in biological research}. New York, W.H. Freeman and Company.
#' 
#' Nakagawa, S. and Schielzeth, H. (2010) \emph{Repeatability for Gaussian and non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
#' 
#' @author Holger Schielzeth  (holger.schielzeth@@ebc.uu.se) & 
#'      Shinichi Nakagawa (shinichi.nakagawa@@otago.ac.nz)
#'      
#' @seealso \link{rpt.aov}, \link{rpt.remlLMM}, \link{rpt.mcmcLMM}, \link{rpt}, \link{print.rpt} 
#' 
#' @examples  
#' # repeatability for male breeding success on a transformed scale
#'   data(Fledglings)
#'   Fledglings$sqrtFledge <- sqrt(Fledglings$Fledge)
#'   attach(Fledglings)
#'   (rpt.Fledge <- rpt.corr(sqrtFledge, MaleID, nboot=10, npermut=10))  # reduced number of iterations
#'   detach(Fledglings)
#' @keywords models
#' 

rpt.corr <- function(y, groups, CI=0.95, nboot=1000, npermut=1000) {
	# initial checks
	if(length(y)!= length(groups)) 
		stop("y and group have to be of equal length")
	if(nboot < 0)	nboot <- 0
	if(npermut < 1) npermut <- 1
	if(any(is.na(y))) 
		stop("missing values in y ")
	if(any(is.na(groups)))
		stop("missing values in groups ")
	if(!all(table(groups)==2))     
		stop("not exactly two data points per group")
	# preparation
	sortix <- sort.int(as.numeric(groups),index.return=TRUE)$ix
	y1     <- y[sortix][seq(1,length(y),by=2)]
	y2     <- y[sortix][seq(2,length(y),by=2)]
	k      <- length(unique(groups))
	# functions: point estimates of R
	R.pe  <- function(y1, y2, k) {
		y <- c(y1, y2)
		R <- 1/(k-1)*sum((y1-mean(y))*(y2-mean(y))) / var(y)
		return(R) 
	}	
	# point estimation according to equations 4 and 5
	R      <- R.pe(y1, y2, k)
	# confidence interval estimation according to equation 6 and 7
	bootstr <- function(y1, y2, k) {
		samp<- sample(1:k, k, replace=TRUE)
		R.pe(y1[samp], y2[samp], k)
	}
	if(nboot > 0)
		R.boot <- replicate(nboot, bootstr(y1, y2, k), simplify=TRUE) 
	else
		R.boot <- NA
	CI.R   <- quantile(R.boot, c((1-CI)/2, 1-(1-CI)/2), na.rm=TRUE)
	se     <- sd(R.boot, na.rm=TRUE)
	# significance test by permutation
	permut   <- function(y1, y2, k) {
		samp <- sample(1:k, k)
		R.pe(y1[samp], y2, k)
	}
	if(npermut > 1) {
		R.permut <- c(R, replicate(npermut-1, permut(y1, y2, k), simplify=TRUE))
		P.permut <- sum(R.permut >= R)/npermut
	}
	else {
		R.permut <- R
		P.permut <- NA
	}
	# return of results
	res <- list(call=match.call(), datatype="Gaussian", method="corr", CI=CI, 
				R=R, se=se, CI.R=CI.R, P=P.permut, R.permut=R.permut) 
	class(res) <- "rpt"
	return(res) 
}