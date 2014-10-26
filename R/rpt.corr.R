rpt.corr <- function(y, groups, CI=0.95, nboot=1000, npermut=1000) {
	# initial checks
	if(length(y)!= length(groups)) 
		stop("y and group have to be of equal length")
	if(nboot < 0)	nboot <- 0
	if(npermut < 1) npermut <- 1
	if(any(is.na(y))) 
		stop("missing values in y ")
	if(any(is.na(groups)))
		stop("missing values in groups ")
	if(!all(table(groups)==2))     
		stop("not exactly two data points per group")
	# preparation
	sortix <- sort.int(as.numeric(groups),index.return=TRUE)$ix
	y1     <- y[sortix][seq(1,length(y),by=2)]
	y2     <- y[sortix][seq(2,length(y),by=2)]
	k      <- length(unique(groups))
	# functions: point estimates of R
	R.pe  <- function(y1, y2, k) {
		y <- c(y1, y2)
		R <- 1/(k-1)*sum((y1-mean(y))*(y2-mean(y))) / var(y)
		return(R) 
	}	
	# point estimation according to equations 4 and 5
	R      <- R.pe(y1, y2, k)
	# confidence interval estimation according to equation 6 and 7
	bootstr <- function(y1, y2, k) {
		samp<- sample(1:k, k, replace=TRUE)
		R.pe(y1[samp], y2[samp], k)
	}
	if(nboot > 0)
		R.boot <- replicate(nboot, bootstr(y1, y2, k), simplify=TRUE) 
	else
		R.boot <- NA
	CI.R   <- quantile(R.boot, c((1-CI)/2, 1-(1-CI)/2), na.rm=TRUE)
	se     <- sd(R.boot, na.rm=TRUE)
	# significance test by permutation
	permut   <- function(y1, y2, k) {
		samp <- sample(1:k, k)
		R.pe(y1[samp], y2, k)
	}
	if(npermut > 1) {
		R.permut <- c(R, replicate(npermut-1, permut(y1, y2, k), simplify=TRUE))
		P.permut <- sum(R.permut >= R)/npermut
	}
	else {
		R.permut <- R
		P.permut <- NA
	}
	# return of results
	res <- list(call=match.call(), datatype="Gaussian", method="corr", CI=CI, 
				R=R, se=se, CI.R=CI.R, P=P.permut, R.permut=R.permut) 
	class(res) <- "rpt"
	return(res) 
}