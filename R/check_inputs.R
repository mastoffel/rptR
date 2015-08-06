#' Checks for data argument and uses non-standard evaluation
#' 
#' Checks inputs for all 
#'
#'@param y response variable
#'@param groups grouping variable
#'@param data data.frame containing variables
#'
#'@return returns a list containing the response and the grouping variable
#'
check_inputs <- function(y, groups, data = NULL) {
    
    if (is.null(data)) {
            if (is.character(y) | is.character(groups)) {
                    stop("Provide a data argument or vector names (non-character)")
            }
            return(out <- list(y = y, groups = groups))
    }
        
    if (!is.null(data)) {
        if (!is.name(y) & !is.character(y)) y <- substitute(y)
        if (!is.name(groups) & !is.character(groups)) groups <- substitute(groups)
        
        if (is.name(y)) y <- eval(y, data, parent.frame())
        if (is.name(groups)) groups <- eval(groups, data, parent.frame())
        
        if (is.character(y)) y <- data[, y]
        if (is.character(groups)) groups <- data[, groups]
        
    }
    out <- list(y = y, groups = groups)
    
} 
