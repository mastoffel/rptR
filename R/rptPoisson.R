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
#' #' # repeatability for female clutch size over two years.
#' data(BroodParasitism)
#' (rpt.Host <- rptPoisson(formula = OwnClutches ~ (1|FemaleID) + (1|Season), grname =c('FemaleID', "Season"), 
#'                         data = BroodParasitism, nboot=10, npermut=10))
#'                         
#' # reduced number of nboot and npermut iterations
#' 
#'        rpt.remlLMM.adj(Tarsus ~ 1 + (1|Sex) + (1|BirdID), c('Sex', 'BirdID'), 
#'  data=BodySize, nboot=10, npermut=10))
#'  
#' # repeatability for male fledgling success
#' data(Fledglings)
#' (rpt.Fledge <- rpt.poisGLMM.multi(data = Fledglings, Fledge, MaleID, nboot=10, npermut=10))  
#' # reduced number of nboot and npermut iterations
#' 
#' 
#' 
#' nind = 50
#' nrep = 50
#' latmu = 0
#' latbv = 0.3
#' latgv = 0.1
#' latrv = 0.2
#' indid = factor(rep(1:nind, each=nrep))
#' groid = factor(rep(1:nrep, nind))
#' obsid = factor(rep(1:I(nind*nrep)))
#' latim = rep(rnorm(nind, 0, sqrt(latbv)), each=nrep)
#' latgm = rep(rnorm(nrep, 0, sqrt(latgv)), nind)
#' latvals = latmu + latim + latgm + rnorm(nind*nrep, 0, sqrt(latrv))
#' expvals = exp(latvals)
#' obsvals = rpois(nind*nrep, expvals)
#' beta0 = latmu
#' beta0 = log(mean(obsvals))
#' md = data.frame(obsvals, indid, obsid, groid)
#'
#' R_est <- rptPoisson(formula = obsvals ~ (1|indid) + (1|groid), grname = c("indid", "groid"), 
#'                     data = md, nboot = 10, link = "log", npermut = 10, parallel = FALSE)
#' R_est2 <- rptPoisson(formula = obsvals ~ (1|indid), grname = "indid", 
#'                     data = md, nboot = 10, link = "log", npermut = 10, parallel = FALSE)
#' @export
#' 

rptPoisson <- function(formula, grname, data, link = c("log", "sqrt"), CI = 0.95, nboot = 1000, 
        npermut = 1000, parallel = FALSE, ncores = NULL) {
        
        # link
        if (length(link) > 1) link <- link[1]
        if (!(link %in% c("log", "sqrt"))) stop("Link function has to be 'log' or 'sqrt'")
        # observational level random effect
        obsid <- factor(1:nrow(data))
        formula <- update(formula,  ~ . + (1|obsid))
        
        mod <- lme4::glmer(formula, data = data, family = poisson(link = link))
        if (nboot < 0) nboot <- 0
        if (npermut < 1) npermut <- 1

        # point estimates of R
        R_pe <- function(formula, data, grname, peYN = FALSE) {
                
                mod <- lme4::glmer(formula = formula, data = data, family = poisson(link = link))
                
                VarComps <- lme4::VarCorr(mod)
                # find groups
                row_group <- which(as.data.frame(VarComps)[["grp"]] %in% grname)
                # 
                var_a <- as.data.frame(VarComps)[["vcov"]][row_group]
                names(var_a) <- as.data.frame(VarComps)[["grp"]][row_group]
                
                var_e = as.numeric(VarComps$obsid)^2
                
                # intercept on link scale
                beta0 <- unname(lme4::fixef(mod)[1])
                
                # varComps shouldn´t contain obsind here
                VarCompsDf <- as.data.frame(VarComps)
                VarCompsGr <- VarCompsDf[which(VarCompsDf[["grp"]] %in% grname), ]
                
                 if (peYN & any(VarCompsGr$vcov == 0)) {
                         if (nboot > 0){
                         assign("nboot", 0, envir = e1)
                         warning("(One of) the point estimate(s) for the repeatability was exactly 
                                 zero; parametric bootstrapping has been skipped.")
                         }
                 }
        
                if (link == "sqrt") {
                        R_link <- var_a/(var_a + var_e * 0.25)
                        R_org <- NA
                }
                
                if (link = "log") {
                        estdv = log(1/exp(beta0)+1)
                        R_link = var_a /(var_a + var_e +  estdv)
                        EY <- exp(beta0 + (var_e + var_a)/2)
                        R_org <- EY * (exp(var_a) - 1)/(EY * (exp(var_e + var_a) - 1) + 1)
                }
                # check whether that works for any number of var
                R <- as.data.frame(rbind(R_org, R_link))
                return(R)
        }
        
        R <- R_pe(formula, data, grname, peYN = TRUE)
        
        # confidence interval estimation by parametric bootstrapping
        Ysim <- as.matrix(simulate(mod, nsim = nboot))
        bootstr <- function(y, mod, formula, data, grname) {
                data[, names(model.frame(mod))[1]] <- as.vector(y)
                R_pe(formula, data, grname)
        }
        if (nboot > 0 & parallel == TRUE) {
                if (is.null(ncores)) {
                        ncores <- parallel::detectCores() - 1
                        warning("No core number specified: detectCores() is used to detect the number of \n cores on the local machine")
                }
                # start cluster
                cl <- parallel::makeCluster(ncores)
                parallel::clusterExport(cl, "R_pe")
                R_boot <- unname(parallel::parApply(cl, Ysim, 2, bootstr, mod = mod, formula = formula, 
                        data = data, grname = grname))
                parallel::stopCluster(cl)
        }
        if (nboot > 0 & parallel == FALSE) {
                R_boot <- unname(apply(Ysim, 2, bootstr, mod = mod, formula = formula, data = data, 
                        grname = grname))
        }
        if (nboot == 0) {
                # R_boot <- matrix(rep(NA, length(grname)), nrow = length(grname))
                R_boot <- NA
        }
        
        # transform bootstrapping repeatabilities into vectors
        boot_org <- list()
        boot_link <- list()
        for (i in 1:length(grname)) {
                boot_org[[i]] <- unlist(lapply(R_boot, function(x) x["R_org", grname[i]]))
                boot_link[[i]] <- unlist(lapply(R_boot, function(x) x["R_link", grname[i]]))
        }
        names(boot_org) <- grname
        names(boot_link) <- grname
        
        calc_CI <- function(x) {
                out <- quantile(x, c((1 - CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
        }
        
        # CI into data.frame and transpose to have grname in rows
        CI_org <- as.data.frame(t(as.data.frame(lapply(boot_org, calc_CI))))
        CI_link <- as.data.frame(t(as.data.frame(lapply(boot_link, calc_CI))))
        
        # se
        se_org <- as.data.frame(t(as.data.frame(lapply(boot_org, sd))))
        se_link <- as.data.frame(t(as.data.frame(lapply(boot_link, sd))))
        names(se_org) <- "se_org"
        names(se_link) <- "se_link"
        

        # significance test by permutation of residuals
        P_permut <- rep(NA, length(grname))
        
        # significance test by likelihood-ratio-test
        terms <- attr(terms(formula), "term.labels")
        randterms <- terms[which(regexpr(" | ", terms, perl = TRUE) > 0)]
        
        # no permutation test
        if (npermut == 1) {
                R_permut <- R
                P_permut <- NA
        }
        
        # significance test by permutation of residuals
        # nperm argument just used for parallisation
        permut <- function(nperm, formula, mod, dep_var, grname, data) {
                # for binom it will be logit 
                y_perm <- rpois(nrow(data), exp(log(fitted(mod)) + sample(resid(mod))))
                data_perm <- data
                data_perm[dep_var] <- y_perm
                out <- R_pe(formula, data_perm, grname)
                out
        }
        # response variable
        dep_var <- as.character(formula)[2]

        # R_permut <- matrix(rep(NA, length(grname) * npermut), nrow = length(grname))
        P_permut <- data.frame(matrix(NA, nrow = 2, ncol = length(grname)),
                row.names = c("P_permut_org", "P_permut_link")) 
        
        if(parallel == TRUE) {
                if (is.null(ncores)) {
                        ncores <- parallel::detectCores()
                        warning("No core number specified: detectCores() is used to detect the number of \n cores on the local machine")
                }
                # start cluster
                cl <- parallel::makeCluster(ncores)
                parallel::clusterExport(cl, "R_pe")
                R_permut <- parallel::parLapply(cl, 1:(npermut-1), permut, formula=formula, 
                        mod=mod, dep_var=dep_var, grname=grname, data = data)
                parallel::stopCluster(cl)
                
        } else if (parallel == FALSE) {
                R_permut <- lapply(1:(npermut - 1), permut, formula, mod, dep_var, grname, data)
                
        }
        
        # adding empirical rpt 
        R_permut <- c(list(R), R_permut)
        
        # equal to boot
        permut_org <- list()
        permut_link <- list()
        for (i in 1:length(grname)) {
                permut_org[[i]] <- unlist(lapply(R_permut, function(x) x["R_org", grname[i]]))
                permut_link[[i]] <- unlist(lapply(R_permut, function(x) x["R_link", grname[i]]))
        }
        names(permut_org) <- grname
        names(permut_link) <- grname
        
        
#         # reshaping and calculating P_permut
#         R_permut_org <- lapply(R_permut, function(x) x["R_org",])
#         R_permut_link <- lapply(R_permut, function(x) x["R_link",])
#         R_permut_org <- do.call(rbind, R_permut_org)
#         R_permut_link <- do.call(rbind, R_permut_link)
        
        P_permut["P_permut_org", ] <- unlist(lapply(permut_org, function(x) sum(x >= x[1])))/npermut
        P_permut["P_permut_link", ] <- unlist(lapply(permut_link, function(x) sum(x >= x[1])))/npermut
        names(P_permut) <- names(permut_link)
        P_permut
        
                
        ## likelihood-ratio-test
        LRT_mod <- as.numeric(logLik(mod))
        LRT_df <- 1
#         if (length(randterms) == 1) {
#                 formula_red <- update(formula, eval(paste(". ~ . ", paste("- (", randterms, ")"))))
#                 LRT.red <- as.numeric(logLik(lm(formula_red, data = data)))
#                 LRT.D <- as.numeric(-2 * (LRT.red - LRT.mod))
#                 LRT.P <- ifelse(LRT.D <= 0, LRT.df, pchisq(LRT.D, 1, lower.tail = FALSE)/2)
#                 # LR <- as.numeric(-2*(logLik(lm(update(formula, eval(paste('. ~ . ', paste('- (',
#                 # randterms, ')') ))), data=data))-logLik(mod))) P.LRT <- ifelse(LR<=0, 1,
#                 # pchisq(LR,1,lower.tail=FALSE)/2)
#         }
        
        
        for (i in c("LRT_P", "LRT_D", "LRT_red")) assign(i, rep(NA, length(grname)))
        
        for (i in 1:length(grname)) {
                formula_red <- update(formula, eval(paste(". ~ . ", paste("- (1 | ", grname[i], 
                        ")"))))
                LRT_red[i] <- as.numeric(logLik(lme4::glmer(formula = formula_red, data = data, 
                        family = poisson(link = link))))
                LRT_D[i] <- as.numeric(-2 * (LRT_red[i] - LRT_mod))
                LRT_P[i] <- ifelse(LRT_D[i] <= 0, 1, pchisq(LRT_D[i], 1, lower.tail = FALSE)/2)
                # LR <- as.numeric(-2*(logLik(lme4::lmer(update(formula, eval(paste('. ~ . ',
                # paste('- (1 | ', grname[i], ')') ))), data=data))-logLik(mod))) P.LRT[i] <-
                # ifelse(LR<=0, 1, pchisq(LR,1,lower.tail=FALSE)/2)
        }
  
        P <- cbind(LRT_P, t(P_permut))
        
        #Function to calculate a point estimate of overdispersion from a mixed model object
        # from Harrison (2014): Using observation-level random effects to
        # model overdispersion in count data in ecology and evolution, PeerJ
        #     od.point<-function(modelobject){
        #             x<-sum(resid(modelobject,type="pearson")^2)
        #             rdf<-summary(modelobject)$AICtab[5]
        #             return(x/rdf)
        #     }
     
        res <- list(call = match.call(), 
                datatype = "Poisson", 
                link = link,
                CI = CI, 
                R = R, 
                se = cbind(se_org,se_link), 
                CI_emp = list(CI_org = CI_org, CI_link = CI_link), 
                P = P,
                R_boot_link = boot_link, 
                R_boot_org = boot_org,
                R_permut_link = permut_link, 
                R_permut_org = permut_org,
                LRT = list(LRT_mod = LRT_mod, LRT_red = LRT_red, LRT_D = LRT_D, LRT_df = LRT_df, 
                LRT_P = LRT_P), 
                ngroups = unlist(lapply(data[grname], function(x) length(unique(x)))), 
                nobs = nrow(data), mod = mod)
        class(res) <- "rpt"
        return(res)
} 
