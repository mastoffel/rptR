rpt.binomGLMM.add <- function(y, groups, CI=0.95, prior=NULL, verbose=FALSE, ...) {
	# initial checks
	if (is.null(dim(y))) 
		y      <- cbind(y, 1-y)
	if(any(is.na(y))) warning("missing values in y")
	# preparation	
	n          <- rowSums(y)
	groups     <- factor(groups)
	# model fitting
	if(all(n==1)) {
		if(is.null(prior)) prior=list(R=list(V=1,fix=1),G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=25^2)))
		mod    <- MCMCglmm(m ~ 1, random= ~ groups, data=data.frame(m=y[,1],nm=y[,2],groups=groups), prior=prior, family="categorical", verbose=verbose, ...) 
	}
	else {
		if(is.null(prior)) prior=list(R=list(V=1e-10,nu=-1),G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=25^2)))
		mod    <- MCMCglmm(cbind(m, nm) ~ 1, random= ~ groups, data=data.frame(m=y[,1],nm=y[,2],groups=groups), prior=prior, family="multinomial2", verbose=verbose, ...)
	}
	# ezxtraction of posterior distributions
	var.a      <- mod$VCV[,"groups"]
	var.e      <- mod$VCV[,"units"]
	postR.link <- var.a / (var.a + var.e + pi^2 /3)
	beta0      <- mod$Sol
	P          <- exp(beta0) / (1+ exp(beta0))
	postR.org  <- (var.a * P*P / (1+exp(beta0))^2 ) / ((var.a + var.e)* P*P / (1+exp(beta0))^2 +  P*(1-P))
	# point estimate on link and original scale
	R.link     <- posterior.mode( postR.link )
	R.org      <- posterior.mode( postR.org )
	# credibility interval estimation from paterior distribution
	CI.link    <- coda::HPDinterval(postR.link, CI)[1,]
	CI.org     <- coda::HPDinterval(postR.org, CI)[1,]
	se.link    <- sd(postR.link)
	se.org     <- sd(postR.org)
	# P value equivilents that are not appropriate
	P.link     <- NA
	P.org      <- NA
	# return of results
	res        <- list(call=match.call(), datatype="binomial", method="MCMC", CI=CI,
			      R.link=R.link, se.link=se.link, CI.link=CI.link, P.link=P.link, 
				  R.org = R.org, se.org=se.org, CI.org=CI.org, P.org=P.org,
				  R.post=list(R.link=as.vector(postR.link), R.org=as.vector(postR.org)) )
	class(res) <- "rpt"
	return(res) 
}
