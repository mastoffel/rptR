rpt.poisGLMM.multi = function(y, groups, link=c("log", "sqrt"), CI=0.95, nboot=1000, npermut=1000) {
	# initial checks
	if(length(y) != length(groups)) 
		stop("y and group hav to be of equal length")
	if(nboot < 0)   nboot <- 0
	if(npermut < 1) npermut <- 1
	if(length(link) > 1)
		link   <- link[1]
	if(link != "log" &  link != "sqrt") 
		stop("inappropriate link (has to be 'log' or 'sqrt')")
	if(any(is.na(y))) {
		warning("missing values in y are removed")
		groups <- groups[!is.na(y)]
		y      <- y[!is.na(y)]
	}
	# preparation
	groups <- factor(groups)
	N <- length(y)
	k <- length(levels(groups))
	# functions
	pqlglmm.pois.model <- function(y, groups, link, returnR=TRUE) {
		mod     <-  glmmPQL(y ~ 1,random=~1|groups,  family=quasipoisson(link=eval(link)), verbose=FALSE) 
		VarComp <- nlme::VarCorr(mod)
		beta0   <- as.numeric(mod$coefficients$fixed)
		omega   <- (as.numeric(VarComp[2,1]))
		var.a   <- (as.numeric(VarComp[1,1]))
		if (link=="log") {
			R.link  <- var.a/(var.a + omega*log(1/exp(beta0)+1))
			EY 		<- exp(beta0+var.a/2)
			R.org 	<- EY*(exp(var.a)-1)/(EY*(exp(var.a)-1)+omega) 
		}
		if (link=="sqrt") {
			R.link  <- var.a/(var.a + omega*0.25)
			R.org 	<- NA 
		}	
		if(returnR) return(list(R.link=R.link, R.org=R.org))
		else return(list(beta0=beta0, omega=omega, var.a=var.a))
	}
	# point estimation according to model 17 equations 18-20
	R   <- pqlglmm.pois.model(y, groups, link)
	# confidence interval estimation by parametric bootstrapping
	bootstr <- function(y, groups, k, N, beta0, var.a, omega, link) {
		groupMeans <- rnorm(k, 0, sqrt(var.a))
		if(link=="log")  mu <- exp(beta0 + groupMeans[groups])
		if(link=="sqrt") mu <- (beta0 + groupMeans[groups])^2
		if (omega<=1)    y.boot <- rpois(N, mu)
			else         y.boot <- rnbinom(N, size=(mu/(omega-1)), mu=mu)
		pqlglmm.pois.model(y.boot, groups, link) 
	}
	if(nboot > 0) {
		mod.ests <- pqlglmm.pois.model(y, groups, link, returnR=FALSE)
		R.boot   <- replicate(nboot, bootstr(y, groups, k, N, mod.ests$beta0, mod.ests$var.a, mod.ests$omega, link), simplify=TRUE)
		R.boot   <- list(R.link = as.numeric(unlist(R.boot["R.link",])), R.org = as.numeric(unlist(R.boot["R.org",]))) 	
	}
	else {
		R.boot   <- list(R.link = NA, R.org = NA)
	}
	CI.link  <- quantile(R.boot$R.link, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
	CI.org   <- quantile(R.boot$R.org, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
	se.link  <- sd(R.boot$R.link,na.rm=TRUE)
	se.org   <- sd(R.boot$R.org,na.rm=TRUE)
	# significance test by randomization
	permut   <- function(y, groups, N, link) {
		samp <- sample(1:N, N)
		pqlglmm.pois.model(y, groups[samp], link) 
	}
	if(npermut > 1) {
		R.permut <- replicate(npermut-1, permut(y, groups, N, link), simplify=TRUE)
		R.permut <- list(R.link = c(R$R.link, unlist(R.permut["R.link",])), R.org = c(R$R.org, unlist(R.permut["R.org",])))
		P.link   <- sum(R.permut$R.link >= R$R.link) / npermut
		P.org    <- sum(R.permut$R.org >= R$R.org) / npermut
	}
	else {
		R.permut <- R
		P.link   <- NA
		P.org    <- NA
	}
	# return of results
	if(mod.ests$omega < 1) 
		warning("omega < 1, therefore CI limits are unreliable")
	res <- list(call=match.call(), datatype="count", method="PQL", link=link, CI=CI,
				R.link = R$R.link, se.link=se.link, CI.link=CI.link, P.link=P.link,
				R.org  = R$R.org, se.org=se.org, CI.org=CI.org, P.org=P.org, 
				omega=mod.ests$omega,
				R.boot = list(R.link=R.boot$R.link, R.org=R.boot$R.org),
				R.permut = list(R.link=R.permut$R.link, R.org=R.permut$R.org) ) 
	class(res) <- "rpt"
	return(res)		
}			
