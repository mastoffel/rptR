#' LMM-based Repeatability Using REML
#' 
#' Calculates repeatability from a linear mixed-effects models fitted by REML (restricted maximum likelihood).
#' 
#' @param y Vector of a response values.
#' @param groups Vector of group identities.
#' @param CI Width of the confidence interval (defaults to 0.95).
#' @param nboot Number of parametric bootstraps for interval estimation. 
#'        Defaults to 1000. Larger numbers of permutations give a better 
#'        asymtotic CI, but may be very time-consuming.
#' @param npermut Number of permutations used when calculating 
#'        asymptotic \emph{P} values (defaults to 1000).
#' 
#' 
#' @return 
#' Returns an object of class rpt that is a a list with the following elements: 
#' \item{datatype}{Response distribution (here: "Gaussian").}
#' \item{method}{Method used to calculate repeatability (here: "REML").}
#' \item{CI}{Width of the confidence interval.}
#' \item{R}{Point estimate for repeatability.}
#' \item{se}{Approximate standard error (\emph{se}) for repeatability. Note that the distribution might not be symmetrical, in which case the \emph{se} is less informative.}
#' \item{CI.R}{Confidence interval for  repeatability.}
#' \item{P}{Approximate \emph{P} value from a significance test based on permutation.}
#' \item{R.boot}{Parametric bootstrap samples for \emph{R}.}
#' \item{R.permut}{Permutation samples for \emph{R}.}
#'
#' @references 
#' Carrasco, J. L. and Jover, L.  (2003). \emph{Estimating the generalized 
#' concordance correlation coefficient through variance components}. Biometrics 59: 849-858.
#'
#' Faraway, J. J. (2006). \emph{Extending the linear model with R}. Boca Raton, FL, Chapman & Hall/CRC.
#' 
#' Nakagawa, S. and Schielzeth, H. (2010) \emph{Repeatability for Gaussian and 
#'              non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
#' 
#' @author Holger Schielzeth  (holger.schielzeth@@ebc.uu.se) & 
#'      Shinichi Nakagawa (shinichi.nakagawa@@otago.ac.nz)
#'      
#' @seealso \link{rpt.mcmcLMM}, \link{rpt.aov}, \link{rpt.corr}, \link{print.rpt}, \link{rpt.remlLMM.adj}
#' 
#' @examples  
#' \dontrun{
#' # repeatability estimation for tarsus length - a very high R
#' data(BodySize)
#' attach(BodySize)
#' (rpt.BS <- rpt.remlLMM(Tarsus, BirdID, nboot=10, npermut=10))   
#' # reduced number of nboot and npermut iterations
#' detach(BodySize)
#'
#' # repeatability estimation for weight (body mass) - a lower R than the previous one
#' data(BodySize)
#' attach(BodySize)
#' (rpt.Weight <- rpt.remlLMM(Weight, BirdID, nboot=10, npermut=10)) 
#' # reduced number of nboot and npermut iterations
#' detach(BodySize)
#'       }
#'       
#' @keywords models
#' 
#' @export
#' 
#' @importFrom nlme VarCorr
#' @importFrom nlme lme
#' 
#' 
rpt.remlLMM <- function(y, groups, CI=0.95, nboot=1000, npermut=1000) {
	# initial checks
	if(length(y)!= length(groups)) stop("y and group are of unequal length")
	if(nboot < 0) 	nboot <- 0
	if(npermut < 1) npermut <- 1
	if(any(is.na(y))) {
		warning("missing values in y are removed")
		groups <- groups[!is.na(y)]
		groups <- factor(groups)
		y      <- y[!is.na(y)] 
	}
	# preparation
	groups <- factor(groups)
    k <- length(unique(groups))
	N <- length(y)
	# functions: point estimates of R
	R.pe <- function(y, groups) {
        varComps <- try(nlme::VarCorr(lme(y ~ 1, random = ~1 | groups)), silent=TRUE)
		if(class(varComps)=="try-error") {
			warning("Convergence problems in lme (most likely during bootstrap or permutation). R is set to NA for this iteration")
			R = NA
		}
		else {
			var.a <- as.numeric(varComps[1, 1])
			var.e <- as.numeric(varComps[2, 1])
			R <- var.a/(var.a + var.e)
		}
		return(R) 
	}
	# point estimation according to model 8 and equation 9
	R <- R.pe(y, groups)
	# confidence interval estimation by parametric bootstrapping
	bootstr <- function(y, groups, k, N, beta0, var.a, var.e) {
		y.boot <- beta0 + rnorm(k, 0, sqrt(var.a))[groups] + rnorm(N, 0, sqrt(var.e))
		R.pe(y.boot, groups) 
	}
	if(nboot > 0) {
		mod <- try(lme(y ~ 1, random = ~1 | groups), silent=TRUE)
		if(class(mod)=="try-error") {
			warning("Convergence problems in lme. Model is refitted.")
			mod <- lme(y ~ 1, random = ~1 | groups)                
		}
		beta0    <- as.numeric(summary(mod)$tTable[,1])
		varComps <- nlme::VarCorr(mod)
		var.a    <- as.numeric(varComps[1,1])
		var.e    <- as.numeric(varComps[2,1])	
		R.boot   <- replicate(nboot, bootstr(y, groups, k, N, beta0, var.a, var.e), simplify=TRUE)
	}
	else
		R.boot   <- NA
	CI.R     <- quantile(R.boot, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
	se       <- sd(R.boot, na.rm=TRUE)
	# significance test by likelihood-ratio-test
	LR       <- as.numeric(-2*(logLik(lm(y~1))-logLik(mod)))
	P.LRT    <- ifelse(LR<=0, 1, pchisq(LR,1,lower.tail=FALSE)/2)
	# significance test by permutation
	permut <- function(y, groups, N) {
		samp <- sample(1:N, N)
		R.pe(y, groups[samp]) 
	}
	if(npermut > 1) {
		R.permut <- c(R, replicate(npermut-1, permut(y, groups, N), simplify=TRUE))
		P.permut <- sum(R.permut >= R)/npermut
	}
	else {
		R.permut = R
		P.permut <- NA
	}
	# return of results
	res  <- list(call=match.call(), datatype="Gaussian", method="LMM.REML", CI=CI, 
				 R=R, se=se, CI.R=CI.R, 
				 P = c(P.LRT=P.LRT, P.permut=P.permut), 
				 R.boot=R.boot, R.permut=R.permut )
	class(res) <- "rpt"
	return(res) 
}
