rpt.poisGLMM.add <- function(y, groups, CI=0.95, prior=NULL, verbose=FALSE, ...) {
	# initial checks
	if(any(is.na(y))) warning("missing values in y")
	# preparation
	groups     <- factor(groups)
	# model fitting
	if(is.null(prior)) prior=list(R=list(V=1e-10,nu=-1),G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=25^2)))
	mod        <- MCMCglmm(y ~ 1, random=~groups, family="poisson", data=data.frame(y=y,groups=groups), prior=prior, verbose=verbose, ...)
	# ezxtraction of posterior distributions
	var.a      <- mod$VCV[,"groups"]
	var.e      <- mod$VCV[,"units"]
	beta0      <- mod$Sol
	postR.link <- var.a/(var.a + var.e+log(1/exp(beta0)+1))
	EY 		   <- exp(beta0+(var.e+var.a)/2)
	postR.org  <- EY*(exp(var.a)-1)/(EY*(exp(var.e+var.a)-1)+1)
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
	res 	   <- list(call=match.call(), datatype="count", method="MCMC", CI=CI,
			      R.link=R.link, se.link=se.link, CI.link=CI.link, P.link=P.link, 
				  R.org = R.org, se.org=se.org, CI.org=CI.org, P.org=P.org,
				  R.post=list(R.link=as.vector(postR.link), R.org=as.vector(postR.org)) )
	class(res) <- "rpt"
	return(res) 
}
