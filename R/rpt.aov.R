rpt.aov <- function(y, groups, CI=0.95, npermut=1000) {	
	# initial checks
    if(length(y) != length(groups)) stop("y and group are not of equal length")
	if(npermut < 1) npermut <- 1
	# preparation
	groups <- factor(groups)
    k  <- length(levels(groups))
	N  <- length(y)
    n0 <- 1/(k-1)*(N - sum(table(groups)^2)/N)
	# functions: point estimates of R
	R.pe <- function(y, groups, n0) {
		gm   <- mean(y)
		MSa  <- sum(tapply(y, groups, function(x) (mean(x)-gm)^2 * length(x)))   / (k-1)
		MSw  <- sum(tapply(y, groups, function(x) sum((x - mean(x))^2)))  / (N-k)
		R    <- ((MSa-MSw)/n0 ) / ( (MSa-MSw)/n0 + MSw)
		return(R) 
	}
	# point estimation according to equations 4 and 5
	R        <- R.pe(y, groups, n0)
	# confidence interval estimation according to equation 6 and 7
	se       <- sqrt(( 2*(N-1)*(1-R)^2 * (1 + (n0-1)*R)^2) / (n0^2 * (N-k)*(k-1)))
	CI.R     <- R + c(1,-1)* qt((1-CI)/2,k-1)*se
	# significance test from ANOVA
	P.aov    <- anova(lm(y ~ groups))[5][1,1]
	# significance test by permutation
	permut <- function(y, groups, N, n0) {
		sampy <- sample(y, N)
		return(R.pe(sampy, groups, n0))
	}
	if(npermut > 1) {
		R.permut <- c(R, replicate(npermut-1, permut(y, groups, N, n0), simplify=TRUE))
		P.permut <- sum(R.permut >= R)/npermut
	}
	else {
		R.permut <- R
		P.permut <- NA
	}
	# return of results
	res <- list(call=match.call(), datatype="Gaussian", method="ANOVA", CI=CI, R=R, se=se, CI.R=CI.R, P=c(P.aov=P.aov, P.permut=P.permut), R.permut=R.permut) 
	class(res) <- "rpt"
	return(res) 
}
