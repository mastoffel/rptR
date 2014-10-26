rpt.binomGLMM.multi <- function(y, groups, link=c("logit", "probit"), CI=0.95, nboot=1000, npermut=1000) {
	# initial checks
	if (is.null(dim(y))) 
		y <- cbind(y, 1-y)
	if (nrow(y) != length(groups)) 
		stop("y and group have to be of equal length")
	if(nboot < 0) 	nboot <- 0
	if(npermut < 1) npermut <- 1
	if(length(link) > 1)
		link   <- link[1]
	if(link != "logit" & link != "probit") 
		stop("inappropriate link (has to be 'logit' or 'probit')")
	if(any(is.na(y))) {
		warning("missing values in y are removed")
		groups <- groups[-rowSums(is.na(y)) > 0]
		y      <- y[-rowSums(is.na(y)) > 0,]
	}
	# preparation
	groups <- factor(groups)
	n <- rowSums(y)
	N <- nrow(y)
	k <- length(levels(groups))
	# functions
	pqlglmm.binom.model <- function(y, groups, n, link, returnR=TRUE) {
		if(all(n==1)) mod <- glmmPQL(y ~ 1, random=~1|groups,  family=binomial(link=eval(link)), verbose=FALSE)
		else          mod <- glmmPQL(y ~ 1, random=~1|groups,  family=quasibinomial(link=eval(link)), verbose=FALSE)	
		VarComp  <- nlme::VarCorr(mod)
		beta0    <- as.numeric(mod$coefficients$fixed)
		if(all(n==1)) omega <- 1
			else      omega <- (as.numeric(VarComp[2,1]))
		var.a  <- (as.numeric(VarComp[1,1]))
		if (link=="logit") {
			R.link  <- var.a / (var.a+omega*pi^2 /3)
			P 		<- exp(beta0) / (1+ exp(beta0))
			R.org   <- (var.a * P*P / (1+exp(beta0))^2 ) / (var.a * P*P / (1+exp(beta0))^2 + omega * P*(1-P)) 
		}
		if (link=="probit") {
			R.link  <- var.a / (var.a+omega)
			R.org   <- NA 
		}
		if(returnR) return(list(R.link=R.link, R.org=R.org))
		else return(list(beta0=beta0, omega=omega, var.a=var.a)) 
	}
	# point estimate
	R <- pqlglmm.binom.model(y, groups, n, link)
	# confidence interval by parametric bootstrapping
	bootstr <- function(y, groups, k, N, n, beta0, var.a, omega, link) {
		groupMeans <- rnorm(k, 0, sqrt(var.a))
		p.link     <- beta0 + groupMeans[groups]
		if(link=="logit") p <- exp(p.link)/(1+exp(p.link))
		if(link=="probit") p <- probit(p.link, inverse=TRUE)
		if(all(n==1))
			m <- rbinom(N, 1, p)   # binomial model
		else {
			rho <- (omega-1)/(n-1)
			rho[rho<=0] <-0				# underdispersion is ignored
			rho[rho>=1] <- 9e-10      # strong overdispersion is forced to have rho ~ 1
#			if(rho==0)
#				m <- rbinom(N,n,p)
#			else {
				m <- rbetabinom(N,n,p,rho)    # or p * rho
#			}
		}
		pqlglmm.binom.model(cbind(m, n-m), groups, n, link) 
	}
	if(nboot > 0) {
		mod.ests <- pqlglmm.binom.model(y, groups, n, link, returnR=FALSE)
		R.boot   <- replicate(nboot, bootstr(y, groups, k, N, n, mod.ests$beta0, mod.ests$var.a, mod.ests$omega, link), simplify=TRUE)
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
		pqlglmm.binom.model(y, groups[samp], n, link) 
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
	if(mod.ests$omega<1) warning("omega < 1, therefore CI are unreliable")
	res <- list(call=match.call(), 
				datatype="binomial", method="PQL", link=link, CI=CI, 
				R.link=R$R.link, se.link=se.link, CI.link=CI.link, P.link=P.link,
				R.org=R$R.org, se.org=se.org, CI.org=CI.org, P.org=P.org, 
				omega=mod.ests$omega,
				R.boot = list(R.link=R.boot$R.link, R.org=R.boot$R.org),
				R.permut = list(R.link=R.permut$R.link, R.org=R.permut$R.org) ) 
	class(res) <- "rpt"
	return(res) 
}			
