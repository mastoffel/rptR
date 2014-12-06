#' LMM-based Repeatability Using REML
#' 
#' Calculates repeatability from a linear mixed-effects models fitted by REML (restricted maximum likelihood).
#' 
#' @param formula Formula as used e.g. by \link{lmer}. The grouping factor of
#'        interest needs to be included as a random effect, e.g. '(1|groups)'.
#'        Covariates and additional random effects can be included to estimate adjusted repeatabilities.
#' @param grname A character string or vector of character strings giving the
#'        name(s) of the grouping factor(s), for which the repeatability should
#'        be estimated. Spelling needs to match the random effect names as given in \code{fromula}.
#' @param data A dataframe that contains the variables included in the formula argument.
#' @param CI Width of the confidence interval (defaults to 0.95).
#' @param nboot Number of parametric bootstraps for interval estimation.
#'        Defaults to 1000. Larger numbers of permutations give a better
#'        asymtotic CI, but may be very time-consuming.
#' @param npermut Number of permutations used when calculating 
#'        asymptotic \emph{P} values (defaults to 1000). Currently not in use!
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
#' @seealso \link{rpt.mcmcLMM}, \link{print.rpt}, \link{rpt}, \link{rpt.adj}
#' 
#' @examples  
#' 
#' 
#' # repeatability estimation for tarsus length - a very high R
#' data(BodySize)
#' (rpt.BS <- rpt.remlLMM.adj(Tarsus ~ Sex + (1|BirdID), "BirdID", 
#'  data=BodySize, nboot=10, npermut=10))
#' # reduced number of nboot and npermut iterations

#' # repeatability estimation for weight (body mass) - a lower R than the previous one
#' data(BodySize)
#' (rpt.Weight <- rpt.remlLMM.adj(Weight ~ Sex + (1|BirdID), "BirdID", data=BodySize,
#'  nboot=10, npermut=10))
#' # reduced number of nboot and npermut iterations
#' 
#'       
#' @keywords models
#' 
#' @export
#' 
# @importFrom lme4 lmer
# @importFrom lme4 VarCorr
# @importFrom arm sim
 
rpt.remlLMM.adj = function(formula, grname, data, CI=0.95, nboot=1000, npermut=1000) {
	mod         <- lme4::lmer(formula, data=data)
	if(nboot < 0) 	nboot <- 0
	if(npermut < 1) npermut <- 1
	e1 = environment()
	# point estimates of R
	R.pe <- function(formula, data, grname, peYN=FALSE) {
		mod.fnc = lme4::lmer(formula, data)
		varComps <- lme4::VarCorr(mod.fnc)
		if(peYN & any(varComps==0) & nboot > 0) {
			assign("nboot", 0, envir=e1)
			warning("(One of) the point estimate(s) for the repeatability was exactly zero; parametric bootstrapping has been skipped.")
		}
		var.a    <- as.numeric(varComps[grname])
		var.p    <- sum(as.numeric(varComps)) + attr(varComps, "sc")^2
		#var.e    <- as.numeric(attr(varComps, "sc")^2)
		R        <- var.a / var.p
		return(R) 
	}
	R <- R.pe(formula, data, grname, peYN=TRUE)
	names(R) = grname
	# confidence interval estimation by parametric bootstrapping
        
	Ysim <- as.matrix(simulate(mod, nsim = nboot))
	bootstr <- function(y, mod, formula, data, grname) {
		data[,names(model.frame(mod))[1]] = as.vector(y)
		R.pe(formula, data, grname)
	}
	if(nboot > 0){ 
                R.boot   <- unname(apply(Ysim, 2, bootstr, mod = mod, formula = formula, data = data, grname = grname))
		} else {
                  R.boot <- matrix(rep(NA, length(grname)), nrow=length(grname))
		}
	if(length(grname) == 1) {
		CI.R     <- quantile(R.boot, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
		se <- sd(R.boot)
		names(se) = grname 
	} else {
		CI.R     <- t(apply(R.boot, 1, function(x) { quantile(x, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)}))	
		se       <- apply(R.boot, 1, sd)
		rownames(R.boot) = grname
		rownames(CI.R) = grname
		names(se) = grname
	}
	# significance test by permutation
	P.permut <- rep(NA, length(grname))
	# significance test by likelihood-ratio-test
	terms = attr(terms(formula), "term.labels")
	randterms = terms[which(regexpr(" | ",terms,perl=TRUE)>0)]
	if(length(randterms)==1) {
		LR       <- as.numeric(-2*(logLik(lm(update(formula, eval(paste(". ~ . ", paste("- (", randterms, ")") ))), data=data))-logLik(mod)))
		P.LRT    <- ifelse(LR<=0, 1, pchisq(LR,1,lower.tail=FALSE)/2)
	}
	if(length(randterms)>1) {
		P.LRT = rep(NA, length(grname))
		for(i in 1:length(grname)) {
			LR       <- as.numeric(-2*(logLik(lmer(update(formula, eval(paste(". ~ . ", paste("- (1 | ", grname[i], ")") ))), data=data))-logLik(mod)))
			P.LRT[i] <- ifelse(LR<=0, 1, pchisq(LR,1,lower.tail=FALSE)/2)
		}
	}
	# preparing results
	P = matrix(c(P.LRT, P.permut),ncol=2,byrow=FALSE)
	colnames(P) = c("P.LRT", "P.permut")
	rownames(P) = grname
	res  = list(datatype="Gaussian", 
				method="LMM.REML", 
				CI=CI,
				R=R, se=se, CI.R=CI.R, 
				P = P,
				R.boot=R.boot, R.permut=NA,
				mod = mod	)
	class(res) <- "rpt"
	return(res)
}
