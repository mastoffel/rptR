#' Checks for data argument and uses non-standard evaluation

check_inputs <- function(y, groups, data = NULL) {
        
        if (is.null(data)) return(out <- list(y = y, groups = groups))
        if (!is.null(data)) {
                y <- eval(substitute(y), data, parent.frame())
                groups <- eval(substitute(groups), data, parent.frame())
                if(is.character(y) & (length(y) <= 2) & is.character(groups) & (length(groups) == 1)) {
                        y <- data[, y]
                        groups <- data[, groups]
                }
        }
        out <- list(y = y, groups = groups)
        
}