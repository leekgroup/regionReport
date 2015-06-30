#' Set an advanced argument
#'
#' @param name Name of the advanced argument to look for in ...
#' @param value The default value of the advanged argument
#' @keywords internal 
.advanced_argument <- function(name, value, ...) {
    args <- list(...)
    if(!name %in% names(args)) {
        return(value)
    } else {
        return(args[[name]])
    }
}


#' temporarily evaluate an expression in a directory
#'
#' Temporarily evaluate an expression in a directory, then set the directory
#' back to the original.
#'
#' @param dir a directory to perform an expression within
#' @param expr expression to evaluate
#'
#' @details See here: http://plantarum.ca/code/setwd-part2/
with_wd <- function(dir, expr) {
    wd <- getwd()
    on.exit(setwd(wd))
    setwd(dir)
    eval(expr, envir = parent.frame())
}
