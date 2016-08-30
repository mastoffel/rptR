#' Captures and suppresses (still to find out why) warnings of an expression
#' 
#' This function is used within rptR to capture lme4 model fitting warnings in the
#' bootstrap and permutation procedures.
#' 
#' @param An expression, such as the sequence of code used by rptR to calculate
#' bootstrap or permutation estimates
#' 
#' @keywords internal 


with_warnings <- function(expr) {
        myWarnings <- NULL
        wHandler <- function(w) {
                myWarnings <<- c(myWarnings, list(w))
                invokeRestart("muffleWarning")
        }
        val <- withCallingHandlers(expr, warning = wHandler)
        list(warnings = myWarnings)
} 


