#' Calculate overdispersion
#'
#' @param model merMod object
#' @return dispersion parameter
#' @details From http://glmm.wikidot.com/faq

overdisp <- function(model) {
                ## number of variance parameters in 
                ##   an n-by-n variance-covariance matrix
                vpars <- function(m) {
                        nrow(m)*(nrow(m)+1)/2
                }
                model_df <- sum(sapply(lme4::VarCorr(model),vpars))+length(lme4::fixef(model))
                rdf <- nrow(model.frame(model))-model_df
                rp <- stats::residuals(model,type="pearson")
                pearson_chisq <- sum(rp^2)
                overdisp <- pearson_chisq/rdf
}

