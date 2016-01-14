#' LMM-based Repeatability Using REML
#' 
#' Calculates repeatability from a linear mixed-effects models fitted by REML (restricted maximum likelihood).
#' 
#' @param formula Formula as used e.g. by \link{glmer}. The grouping factor of
#'        interest needs to be included as a random effect, e.g. '(1|groups)'.
#'        Covariates and additional random effects can be included to estimate adjusted repeatabilities.
#' @param grname A character string or vector of character strings giving the
#'        name(s) of the grouping factor(s), for which the repeatability should
#'        be estimated. Spelling needs to match the random effect names as given in \code{formula}.
#' @param data A dataframe that contains the variables included in the formula argument.
#' @param CI Width of the confidence interval (defaults to 0.95).
#' @param nboot Number of parametric bootstraps for interval estimation.
#'        Defaults to 1000. Larger numbers of permutations give a better
#'        asymtotic CI, but may be very time-consuming.
#' @param npermut Number of permutations used when calculating 
#'        asymptotic \emph{P} values (defaults to 1000). Currently not in use!
#' @param parallel If TRUE, bootstraps will be distributed. 
#' @param ncores Specify number of cores to use for parallelization. On default,
#'        all cores but one are used.
#' 
#' @return 
#' Returns an object of class rpt that is a a list with the following elements: 
#' \item{datatype}{Response distribution (here: 'Gaussian').}
#' \item{method}{Method used to calculate repeatability (here: 'REML').}
#' \item{CI}{Width of the confidence interval.}
#' \item{R}{Point estimate for repeatability.}
#' \item{se}{Approximate standard error (\emph{se}) for repeatability. Note that the distribution might not be symmetrical, in which case the \emph{se} is less informative.}
#' \item{CI.R}{Confidence interval for  repeatability.}
#' \item{P}{Approximate \emph{P} value from a significance test based on likelihood-ratio.}
#' \item{P.permut}{Approximate \emph{P} value from a significance test based on permutation of residuals.}
#' \item{R.boot}{Parametric bootstrap samples for \emph{R}.}
#' \item{R.permut}{Permutation samples for \emph{R}.}
#' \item{LRT}{Vector of Likelihood-ratios for the model(s) and the reduced model(s), and \emph{P} value(s) and degrees of freedom for the Likelihood-ratio test} 
#' \item{ngroups}{Number of groups.}
#' \item{nobs}{Number of observations.}
#' \item{mod}{Fitted model.}
#'
#' @references 
#' Carrasco, J. L. and Jover, L.  (2003). \emph{Estimating the generalized 
#' concordance correlation coefficient through variance components}. Biometrics 59: 849-858.
#'
#' Faraway, J. J. (2006). \emph{Extending the linear model with R}. Boca Raton, FL, Chapman & Hall/CRC.
#' 
#' Nakagawa, S. and Schielzeth, H. (2010) \emph{Repeatability for Gaussian and 
#' non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
#' 
#' @author Holger Schielzeth  (holger.schielzeth@@ebc.uu.se),
#'         Shinichi Nakagawa (shinichi.nakagawa@@otago.ac.nz) &
#'         Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'      
#' @seealso \link{rpt}
#' 
#' @examples  
#'  # repeatability estimations for egg dumping (binary data)
#' data(BroodParasitism)
#'                                      
#' (rpt.BroodPar <- rptBinary(cbpYN ~ 1 + (1|FemaleID), "FemaleID", data = BroodParasitism, 
#'                                      nboot=10, npermut=10))
#' 
#' 
#' @export
#' 

rptBinary <- function(formula, grname, data, link = c("logit", "probit"), CI = 0.95, nboot = 1000, npermut = 1000, 
        parallel = FALSE, ncores = NULL) {
        
        
        # link function
        
        if (length(link) > 1) link <- link[1]
        if (link != "logit" & link != "probit") 
                stop("inappropriate link (has to be 'logit' or 'probit')")
        
        mod <- lme4::glmer(formula, data = data, family = binomial(link = eval(link)))
        
        if (nboot < 0) nboot <- 0
        if (npermut < 1) npermut <- 1
    
        # initial checks
        # if (is.null(dim(y))) y <- cbind(y, 1 - y)
        # if (nrow(y) != length(groups)) 
        #        stop("y and group have to be of equal length")
        # if (nboot < 0) nboot <- 0
        # if (npermut < 1) npermut <- 1
        
       
        # preparation
#         groups <- factor(groups)
#         n <- rowSums(y)
#         N <- nrow(y)
#         k <- length(levels(groups))

        
        R_pe <- function(formula, data, grname, link, peYN = FALSE) {
                
                mod <- lme4::glmer(formula, data, family = binomial(link = eval(link)))
                VarComps <- lme4::VarCorr(mod)
                
                if (peYN & any(varComps == 0) & nboot > 0) {
                        assign("nboot", 0, envir = e1)
                        warning("(One of) the point estimate(s) for the repeatability was exactly zero; parametric bootstrapping has been skipped.")
                }
                # intercept on link scale
                beta0 <- unname(fixef(mod)[1])
                # unclear if correct. Dispersion parameter? 
                omega <- sigma(mod)^2
                # find groups
                row_group <- which(as.data.frame(VarComps)[["grp"]] %in% grname)
                # 
                var_a <- as.data.frame(VarComps)[["vcov"]][row_group]
                
                # dispersion parameter (source: https://stat.ethz.ch/pipermail/r-sig-mixed-models/2011q1/015393.html)
                if (link == "logit") {
                        R_link <- var_a/(var_a + omega * pi^2/3)
                        P <- exp(beta0)/(1 + exp(beta0))
                        R_org <- (var_a * P * P/(1 + exp(beta0))^2)/(var_a * P * P/(1 + exp(beta0))^2 + 
                                        omega * P * (1 - P))
                }
                if (link == "probit") {
                        R_link <- var_a/(var_a + omega)
                        R_org <- NA
                }
                
                # name results
                if (length(grname) > 1){
                        names(R_org) <- grname
                        names(R_link) <- grname
                } 
                
                out <- list(R_link = R_link, R_org = R_org)
                # out <- data.frame(R_link = R_link, R_org = R_org)
        }
        

        R <- R_pe(formula, data, grname, link, peYN = TRUE)
        # names(R) <- c(paste(grname, "_link", sep = ""), paste(grname, "_org", sep = ""))
        
        # confidence interval estimation by parametric bootstrapping
        Ysim <- as.matrix(stats::simulate(mod, nsim = nboot))
        bootstr <- function(y, mod, formula, data, grname, link) {
                data[, names(model.frame(mod))[1]] <- as.vector(y)
                R_pe(formula, data, grname, link)
        }
         
        # to do: sort R_link and R_org 
        if (nboot > 0 & parallel == TRUE) {
                if (is.null(ncores)) {
                        ncores <- parallel::detectCores() - 1
                        warning("No core number specified: detectCores() is used to detect the number of \n cores on the local machine")
                }
                # start cluster
                cl <- parallel::makeCluster(ncores)
                R.boot <- unname(parallel::parApply(cl, Ysim, 2, bootstr, mod = mod, formula = formula, 
                        data = data, grname = grname))
                parallel::stopCluster(cl)
        }
        if (nboot > 0 & parallel == FALSE) {
                R.boot <- unname(apply(Ysim, 2, bootstr, mod = mod, formula = formula, data = data, 
                        grname = grname, link = link))
        }
        if (nboot == 0) {
                R.boot <- matrix(rep(NA, length(grname)), nrow = length(grname))
        }
        if (length(grname) == 1) {
                CI.R <- quantile(R.boot, c((1 - CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
                se <- sd(R.boot)
                names(se) <- grname
        } else {
                CI.R <- t(apply(R.boot, 1, function(x) {
                        quantile(x, c((1 - CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
                }))
                se <- apply(R.boot, 1, sd)
                rownames(R.boot) <- grname
                rownames(CI.R) <- grname
                names(se) <- grname
        }
        # significance test by permutation of residuals
        P.permut <- rep(NA, length(grname))
        
        # significance test by likelihood-ratio-test
        terms <- attr(terms(formula), "term.labels")
        randterms <- terms[which(regexpr(" | ", terms, perl = TRUE) > 0)]
        
        # no permutation test
        if (npermut == 1) {
                R.permut <- R
                P.permut <- NA
        }
        
        # significance test by permutation of residuals
        # nperm argument just used for parallisation
        permut <- function(nperm, formula, mod_red, dep_var, grname, i) {
                y_perm <- fitted(mod_red) + sample(resid(mod_red))
                data_perm <- data
                data_perm[dep_var] <- y_perm
                R.pe <- R.pe(formula, data_perm, grname)[i]
        }
        # response variable
        dep_var <- as.character(formula)[2]
        # one random effect, uses lm()
        if (length(randterms) == 1) {
                
                formula_red <- update(formula, eval(paste(". ~ . ", paste("- (", randterms, ")"))))
                mod_red <- lm(formula_red, data = data)
                # R.permut <- c(R, replicate(npermut-1, permut(formula, groups), simplify=TRUE))
                if (parallel == TRUE){
                        if (is.null(ncores)) {
                                ncores <- parallel::detectCores()
                                warning("No core number specified: detectCores() is used to detect the number of \n cores on the local machine")
                        }
                        # start cluster
                        cl <- parallel::makeCluster(ncores)
                        R.permut <- c(R, parallel::parSapply(cl, npermut-1, permut, formula, mod_red, dep_var, grname))
                        parallel::stopCluster(cl)
                        P.permut <- sum(R.permut >= R)/npermut
                        
                } else if (parallel == FALSE) {
                        R.permut <- c(R, replicate(npermut - 1, permut(formula=formula, mod_red=mod_red, dep_var=dep_var, grname=grname), 
                                simplify = TRUE))
                        P.permut <- sum(R.permut >= R)/npermut
                }
        }
        # multiple random effects, uses glmer()
        if (length(randterms) > 1) {
                R.permut <- matrix(rep(NA, length(grname) * npermut), nrow = length(grname))
                P.permut <- rep(NA, length(grname))
                for (i in 1:length(grname)) {
                        formula_red <- update(formula, eval(paste(". ~ . ", paste("- (1 | ", grname[i], 
                                ")"))))
                        mod_red <- lme4::glmer(formula_red, data = data, family = binomial)
                        
                        if(parallel == TRUE) {
                                if (is.null(ncores)) {
                                        ncores <- parallel::detectCores()
                                        warning("No core number specified: detectCores() is used to detect the number of \n cores on the local machine")
                                }
                                # start cluster
                                cl <- parallel::makeCluster(ncores)
                                R.permut[i, ] <- c(R[i], parallel::parSapply(cl, npermut-1, permut, formula, mod_red, dep_var, grname, i))
                                parallel::stopCluster(cl)
                                P.permut[i] <- sum(R.permut[i, ] >= R[i])/npermut
                        } else if (parallel == FALSE) {
                                R.permut[i, ] <- c(R[i], replicate(npermut - 1, permut(formula=formula, mod_red=mod_red, dep_var=dep_var, grname=grname, i=i), simplify = TRUE))
                                P.permut[i] <- sum(R.permut[i, ] >= R[i])/npermut
                        }
                }
        }
        
        ## likelihood-ratio-test
        LRT.mod <- as.numeric(logLik(mod))
        LRT.df <- 1
        if (length(randterms) == 1) {
                formula_red <- update(formula, eval(paste(". ~ . ", paste("- (", randterms, ")"))))
                LRT.red <- as.numeric(logLik(lm(formula_red, data = data)))
                LRT.D <- as.numeric(-2 * (LRT.red - LRT.mod))
                LRT.P <- ifelse(LRT.D <= 0, LRT.df, pchisq(LRT.D, 1, lower.tail = FALSE)/2)
                # LR <- as.numeric(-2*(logLik(lm(update(formula, eval(paste('. ~ . ', paste('- (',
                # randterms, ')') ))), data=data))-logLik(mod))) P.LRT <- ifelse(LR<=0, 1,
                # pchisq(LR,1,lower.tail=FALSE)/2)
        }
        if (length(randterms) > 1) {
                for (i in c("LRT.P", "LRT.D", "LRT.red")) assign(i, rep(NA, length(grname)))
                for (i in 1:length(grname)) {
                        formula_red <- update(formula, eval(paste(". ~ . ", paste("- (1 | ", grname[i], 
                                ")"))))
                        LRT.red[i] <- as.numeric(logLik(lme4::glmer(formula = formula_red, data = data, family = binomial)))
                        LRT.D[i] <- as.numeric(-2 * (LRT.red[i] - LRT.mod))
                        LRT.P[i] <- ifelse(LRT.D[i] <= 0, 1, pchisq(LRT.D[i], 1, lower.tail = FALSE)/2)
                        # LR <- as.numeric(-2*(logLik(lme4::lmer(update(formula, eval(paste('. ~ . ',
                        # paste('- (1 | ', grname[i], ')') ))), data=data))-logLik(mod))) P.LRT[i] <-
                        # ifelse(LR<=0, 1, pchisq(LR,1,lower.tail=FALSE)/2)
                }
        }
        
        
        # preparing results
        # if more than one random term make matrix
        if (nrow(matrix(c(LRT.P, P.permut), ncol = 2, byrow = FALSE)) > 1) {
                P <- matrix(c(LRT.P, P.permut), ncol = 2, byrow = FALSE)
                colnames(P) <- c("P.LRT", "P.permut")
                rownames(P) <- grname
                # else make vector in congruency with rpt.remlLMM
        } else {
                P = c(P.LRT = LRT.P, P.permut = P.permut)
        }
        
        res <- list(call = match.call(), datatype = "Gaussian", method = "LMM.REML", CI = CI, 
                R = R, se = se, CI.R = CI.R, P = P, P.permut = P.permut, R.boot = R.boot, R.permut = R.permut, 
                LRT = list(LRT.mod = LRT.mod, LRT.red = LRT.red, LRT.D = LRT.D, LRT.df = LRT.df, 
                        LRT.P = LRT.P), ngroups = unlist(lapply(data[grname], function(x) length(unique(x)))), 
                nobs = nrow(data), mod = mod)
        class(res) <- "rpt"
        return(res)
} 
