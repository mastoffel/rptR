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
