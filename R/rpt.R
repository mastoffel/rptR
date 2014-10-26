rpt <- function(y, groups, 
			   datatype=c("Gaussian", "binomial", "proportion", "count"),  
			   method=c("corr", "ANOVA", "REML", "MCMC", "GLMM.add", "GLMM.multi"),  
			   link=c("logit", "probit", "log", "sqrt"),
			   CI=0.95, nboot=1000, npermut=1000) {
	if(datatype=="Gaussian") {
		if(length(method)>1) {
			warning("Linear mixed model fitted by REML used by default. Change using argument 'method', if required ('corr', 'ANOVA', 'REML' and 'MCMC' allowed for Gaussian data).")
			method<-"REML" 
		}
		if (method=="REML")  return(rpt.remlLMM(y, groups, CI=CI, nboot=nboot, npermut=npermut))
		if (method=="MCMC")  return(rpt.mcmcLMM(y, groups, CI=CI))
		if (method=="ANOVA") return(rpt.aov(y, groups, CI=CI, npermut=npermut))	
		if (method=="corr")  return(rpt.corr(y, groups, CI=CI, nboot=nboot, npermut=npermut)) 
	}
	if(datatype=="binomial" | datatype=="proportion") {
		if(length(method)>1) {
			warning("Generalised linear mixed model with multiplicative overdispersion fitted by PQL used by default. Change using argument 'method', if required ('GLMM.add' and 'GLMM.multi' allowed for Binomial data).")
			method<-"GLMM.multi" 
		}
		if (method=="GLMM.multi") return(rpt.binomGLMM.multi(y, groups, link, CI=CI, nboot=nboot, npermut=npermut))
		if (method=="GLMM.add") return(rpt.binomGLMM.add(y, groups, CI=CI))
	}
	if(datatype=="count") {
		if(length(method)>1) {
			warning("Generalised linear mixed model with multiplicative overdispersion fitted by PQL used by default. Change using argument 'method', if required ('GLMM.add' and 'GLMM.multi' allowed for count data).")
			method<-"GLMM.multi"
		}
		if(length(link)>1) {
			link="log"
			warning("Log link will be used.")
		}
		if (method=="GLMM.multi") return(rpt.poisGLMM.multi(y, groups, link, CI=CI, nboot=nboot, npermut=npermut)) 
		if (method=="GLMM.add")   return(rpt.poisGLMM.add(y, groups, CI=CI)) 
	} 
}