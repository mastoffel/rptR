print.rpt <- function(x, ...) {
	if(x$datatype=="Gaussian" & length(x$P)==1 & length(x$R)==1) {
		cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n",
			"R  = ", round(x$R,3), "\n",
			"SE = ", round(x$se,3), "\n",
			"CI = [", round(x$CI.R[1],3), ", ", round(x$CI.R[2],3), "]", "\n",
			"P  = ", signif(x$P, 3), "\n\n", 
			sep="")  		
	} 
	if(x$datatype=="Gaussian" & length(x$P)>1 & length(x$R)==1) {
		cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n",
			"R  = ", round(x$R,3), "\n",
			"SE = ", round(x$se,3), "\n",
			"CI = [", round(x$CI.R[1],3), ", ", round(x$CI.R[2],3), "]", "\n",
			"P  = ", signif(x$P[1], 3), " [", attr(x$P, "names")[1], "]", "\n", 
			"     ", signif(x$P[2], 3), " [", attr(x$P, "names")[2], "]", "\n\n", 
			sep="")
	} 
	if(x$datatype=="Gaussian" & length(x$P)==1 & length(x$R)>1) {
		print("h3")
		cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n")
		for(i in 1: length(x$R)) {
			cat("Repeatability for ", names(x$R)[i], "\n",
			"R  = ", round(x$R[i],3), "\n",
			"SE = ", round(x$se[i],3), "\n",
			"CI = [", round(x$CI.R[i,1],3), ", ", round(x$CI.R[i,2],3), "]", "\n",
			"P  = ", signif(x$P[i], 3), "\n\n", 
			sep="")
		}
	}
	if(x$datatype=="Gaussian" & length(x$P)>1 & length(x$R)>1) {
		cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n")
		for(i in 1: length(x$R)) {
			cat("Repeatability for ", names(x$R)[i], "\n",
			"R  = ", round(x$R[i],3), "\n",
			"SE = ", round(x$se[i],3), "\n",
			"CI = [", round(x$CI.R[i,1],3), ", ", round(x$CI.R[i,2],3), "]", "\n",
			"P  = ", signif(x$P[i,1], 3), " [", attr(x$P, "names")[1], "]", "\n", 
			"     ", signif(x$P[i,2], 3), " [", attr(x$P, "names")[2], "]", "\n\n", 
			sep="")
		}
	}

	if(x$datatype!="Gaussian" & x$method=="PQL") {
		cat("\n", "Repeatability calculation using the ", x$method, " method and ", x$link, "link", "\n\n",
			"Esimated overdistpersion (omega) = ", x$omega, "\n\n",
			"Link scale repeatabilities:","\n",
			"R  = ", round(x$R.link,3), "\n",
			"SE = ", round(x$se.link,3), "\n",
			"CI = [", round(x$CI.link[1],3), ", ", round(x$CI.link[2],3), "]", "\n",
			"P  = ", signif(x$P.link, 3), "\n\n", 
			"Origial scale repeatabilities:","\n",
			"R  = ", round(x$R.org,3), "\n",
			"SE = ", round(x$se.org,3), "\n",
			"CI = [", round(x$CI.org[1],3), ", ", round(x$CI.org[2],3), "]", "\n",
			"P  = ", signif(x$P.org, 3), "\n\n", 
			sep="")  
	}
	if(x$datatype!="Gaussian" & x$method=="MCMC") {
		cat("\n", "Repeatability calculation using the ", x$method, " method", "\n\n",
			"Link scale repeatabilities:","\n",
			"R  = ", round(x$R.link,3), "\n",
			"SE = ", round(x$se.link,3), "\n",
			"CI = [", round(x$CI.link[1],3), ", ", round(x$CI.link[2],3), "]", "\n",
			"P  = ", signif(x$P.link, 3), "\n\n", 
			"Origial scale repeatabilities:","\n",
			"R  = ", round(x$R.org,3), "\n",
			"SE = ", round(x$se.org,3), "\n",
			"CI = [", round(x$CI.org[1],3), ", ", round(x$CI.org[2],3), "]", "\n",
			"P  = ", signif(x$P.org, 3), "\n\n", 
			sep="")  
	}	
}	

