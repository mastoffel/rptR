% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rptGaussian.R
\name{rptGaussian}
\alias{rptGaussian}
\title{LMM-based Repeatability Estimation for Gaussian Data}
\usage{
rptGaussian(
  formula,
  grname,
  data,
  CI = 0.95,
  nboot = 1000,
  npermut = 0,
  parallel = FALSE,
  ncores = NULL,
  ratio = TRUE,
  adjusted = TRUE,
  rptObj = NULL,
  update = FALSE,
  ...
)
}
\arguments{
\item{formula}{Formula as used e.g. by \link[lme4]{lmer}. The grouping factor(s) of
interest needs to be included as a random effect, e.g. '(1|groups)'.
Covariates and additional random effects can be included to estimate adjusted 
repeatabilities.}

\item{grname}{A character string or vector of character strings giving the
name(s) of the grouping factor(s), for which the repeatability should
be estimated. Spelling needs to match the random effect names as given in \code{formula} 
and terms have to be set in quotation marks. The reseved terms "Residual", 
"Overdispersion" and "Fixed" allow the estimation of overdispersion variance, residual 
variance and variance explained by fixed effects, respectively.}

\item{data}{A dataframe that contains the variables included in the \code{formula}
and \code{grname} arguments.}

\item{CI}{Width of the required confidence interval between 0 and 1 (defaults to 
0.95).}

\item{nboot}{Number of parametric bootstraps for interval estimation 
(defaults to 1000). Larger numbers of bootstraps give a better
asymtotic CI, but may be time-consuming. Bootstrapping can be switch off by setting 
\code{nboot = 0}. See also \strong{Details} below.}

\item{npermut}{Number of permutations used when calculating asymptotic p-values 
(defaults to 0). Larger numbers of permutations give a better
asymtotic p-values, but may be time-consuming (in particular when multiple grouping factors
are specified). Permutaton tests can be switch off by setting \code{npermut = 0}. 
See also \strong{Details} below.}

\item{parallel}{Boolean to express if parallel computing should be applied (defaults to FALSE). 
If TRUE, bootstraps and permutations will be distributed across multiple cores.}

\item{ncores}{Specifying the number of cores to use for parallelization. On default,
all but one of the available cores are used.}

\item{ratio}{Boolean to express if variances or ratios of variance should be estimated. 
If FALSE, the variance(s) are returned without forming ratios. If TRUE (the default) ratios 
of variances (i.e. repeatabilities) are estimated.}

\item{adjusted}{Boolean to express if adjusted or unadjusted repeatabilities should be estimated. 
If TRUE (the default), the variances explained by fixed effects (if any) will not
be part of the denominator, i.e. repeatabilities are calculated after controlling for 
variation due to covariates. If FALSE, the varianced explained by fixed effects (if any) will
be added to the denominator.}

\item{rptObj}{The output of a rptR function. Can be specified in combination with update = TRUE
to update bootstraps and permutations}

\item{update}{If TRUE, the rpt object to be updated has to be inputted with the rptObj argument.
The function just updates the permutations and bootstraps, so make sure to specify all other
arguments excactly like for the rpt object specified in rptObj.}

\item{...}{Other parameters for the lmer or glmer call, such as optimizers.}
}
\value{
Returns an object of class \code{rpt} that is a a list with the following elements: 
\item{call}{Function call.}
\item{datatype}{Response distribution (here: 'Gaussian').}
\item{CI}{Coverage of the confidence interval as specified by the \code{CI} argument.}
\item{R}{\code{data.frame} with point estimates for repeatabilities for each grouping factor
      of interest.}
\item{se}{\code{data.frame} with approximate standard errors (\emph{se}) for repeatabilities. 
     Rows repsresent grouping factors of interest. Note that the distribution might not be symmetrical, 
     in which case the \emph{se} is less informative.}
\item{CI_emp}{\code{data.frame} containing the (empirical) confidence intervals for the repeatabilities 
     estiamted based parametric bootstrapping. Each row represents a grouping factor of interest.}
\item{P}{\code{data.frame} with p-values based on likelihood-ratio tests
     (first column) and permutation tests (second column). Each row represents a grouping factor 
     of interest.}
\item{R_boot}{Vector(s) of parametric bootstrap samples for \emph{R}. Each \code{list}
      element respesents a grouping factor.}
\item{R_permut}{Vector(s) of permutation samples for \emph{R}. Each \code{list}
      element represents a grouping factor.}
\item{LRT}{\code{list} with two elements. (1) The likelihood for the full model and a \code{data.frame} 
     called \code{LRT_table} for the reduced model(s), which includes columns
     for the respective grouping factor(s), the likelihood(s) \emph{logL_red}, likelihood ratio(s)
     \emph{LR_D}, p-value(s) \emph{LRT_P} and degrees of freedom \emph{LRT_df}} 
\item{ngroups}{Number of groups for each grouping level.}
\item{nobs}{Number of observations.}
\item{mod}{Fitted model.}
\item{ratio}{Boolean. TRUE, if ratios have been estimated, FALSE, if variances have been estimated}
\item{adjusted}{Boolean. TRUE, if estimates are adjusted}
\item{all_warnings}{\code{list} with two elements. 'warnings_boot' and 'warnings_permut' contain
warnings from the lme4 model fitting of bootstrap and permutation samples, respectively.}
}
\description{
Estimates the repeatability from a general linear mixed-effects models fitted by restricted maximum likelihood (REML).
}
\details{
see details section of \code{\link{rpt}} for details on parametric bootstrapping,
permutation and likelihood-ratio tests.
}
\examples{
 

data(BeetlesBody)

# Note: nboot and npermut are set to 3 for speed reasons. Use larger numbers
# for the real analysis.

# one random effect
rpt_est <- rptGaussian(BodyL ~ (1|Population), grname="Population", 
                   data=BeetlesBody, nboot=3, npermut=3, ratio = FALSE)

# two random effects
rptGaussian(BodyL ~ (1|Container) + (1|Population), grname=c("Container", "Population"), 
                   data=BeetlesBody, nboot=3, npermut=3)
                   
# unadjusted repeatabilities with fixed effects and 
# estimation of the fixed effect variance
rptGaussian(BodyL ~ Sex + Treatment + Habitat + (1|Container) + (1|Population), 
                  grname=c("Container", "Population", "Fixed"), 
                  data=BeetlesBody, nboot=3, npermut=3, adjusted=FALSE)
                  
                  
# two random effects, estimation of variance (instead repeatability)
R_est <- rptGaussian(formula = BodyL ~ (1|Population) + (1|Container), 
            grname= c("Population", "Container", "Residual"),
            data=BeetlesBody, nboot=3, npermut=3, ratio = FALSE)

}
\references{
Carrasco, J. L. & Jover, L.  (2003) \emph{Estimating the generalized 
concordance correlation coefficient through variance components}. Biometrics 59: 849-858.

Nakagawa, S. & Schielzeth, H. (2010) \emph{Repeatability for Gaussian and 
non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
}
\seealso{
\link{rpt}
}
\author{
Holger Schielzeth  (holger.schielzeth@uni-jena.de),
        Shinichi Nakagawa (s.nakagawa@unsw.edu.au) &
        Martin Stoffel (martin.adam.stoffel@gmail.com)
}
