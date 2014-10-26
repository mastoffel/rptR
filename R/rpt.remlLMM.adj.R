rpt.remlLMM.adj = function(formula, grname, data, CI=0.95, nboot=1000, npermut=1000) {
	mod         <- lmer(formula, data=data)
	if(nboot < 0) 	nboot <- 0
	if(npermut < 1) npermut <- 1
	e1 = environment()
	# point estimates of R
	R.pe <- function(formula, data, grname, peYN=FALSE) {
		mod.fnc = lmer(formula, data)
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
	bootstr <- function(mod, formula, data, grname) {
		mod.sim <- arm::sim(mod, n.sim=2)   # for some reason it is not possible to set n.sim=1
		y <- mod@X %*% matrix(mod.sim@fixef[1,]) + t(as.matrix(mod@Zt)) %*% matrix(unlist(mod.sim@ranef)[seq(1, length(unlist(mod.sim@ranef)), by=2)])
		y <- y + rnorm(length(y), 0, attr(lme4::VarCorr(mod), "sc"))
		data[,names(mod@frame)[1]] = as.vector(y)
		R.pe(formula, data, grname)
	}
	if(nboot > 0) R.boot   <- replicate(nboot, bootstr(mod, formula, data, grname), simplify=TRUE)
		else R.boot <- matrix(rep(NA, length(grname)), nrow=length(grname))
	if(length(grname) == 1) {
		CI.R     <- quantile(R.boot, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
		se <- sd(R.boot)
		names(se) = grname }
	else {
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
