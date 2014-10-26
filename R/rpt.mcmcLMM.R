rpt.mcmcLMM <- function(y, groups, CI=0.95, prior=NULL, verbose=FALSE, ...){
	# initial checks
	if(length(y)!= length(groups)) stop("y and group are of unequal length")
	# preparation
	groups <- factor(groups)
	if(is.null(prior)) prior <- list(R=list(V=1,n=10e-2), G=list(G1=list(V=1,n=10e-2)) )
 	# point estimation according to model 8 and equation 9
	mod   <- MCMCglmm(y ~ 1, random=~groups, family="gaussian", data=data.frame(y=y,groups=groups), prior=prior, verbose=verbose, ...)
	var.a <- mod$VCV[,"groups"]
	var.e <- mod$VCV[,"units"]
	postR <- var.a / (var.a + var.e)
	# point estimate
	R     <- posterior.mode( postR )
	# credibility interval estimation from paterior distribution
	CI.R    <- coda::HPDinterval(postR,CI)[1,]
	se 	    <- sd(postR)
	# 'significance test'
	P 	  <- NA
	res = list(call=match.call(), datatype="Gaussian", method="LMM.MCMC", CI=CI, 
				R=R, CI.R=CI.R, se=se, P=P, R.post=postR ) 
	class(res) <- "rpt"
	return(res) 
}
