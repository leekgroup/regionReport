#' Set an advanced argument
#'
#' @param name Name of the advanced argument to look for in ...
#' @param value The default value of the advanged argument
#' @keywords internal
#' @noRd
.advanced_argument <- function(name, value, ...) {
    args <- list(...)
    if (!name %in% names(args)) {
        return(value)
    } else {
        return(args[[name]])
    }
}


#' Temporarily evaluate an expression in a directory
#'
#' Temporarily evaluate an expression in a directory, then set the directory
#' back to the original.
#'
#' @param dir a directory to perform an expression within
#' @param expr expression to evaluate
#'
#' @details See here: http://plantarum.ca/code/setwd-part2/
#' @author Tyler Smith, contributed to regionReport by David Robinson
#' https://github.com/dgrtwo
#' @noRd
with_wd <- function(dir, expr) {
    wd <- getwd()
    on.exit(setwd(wd))
    setwd(dir)
    eval(expr, envir = parent.frame())
}


#' Attempt to load the namespace of one or more packages
#'
#' This function uses requireNamespace to try to load one or more packages.
#' If a package is missing, it will suggest how to install it via Bioconductor
#' before quitting.
#'
#' @param pkg A character vector with the names of the packages to check.
#'
#' @details
#' Updated after feedback from Marcel Ramos at
#' \url{https://github.com/leekgroup/recount/issues/14}
#'
#' @author Leonardo Collado-Torres
#'
#' @noRd
#'
load_check <- function(pkg) {
    ## Check the input
    stopifnot(is.character(pkg))
    ## Try to load the packages
    works <- sapply(pkg, requireNamespace, quietly = TRUE)

    ## If some failed, then give a useful error before quitting
    if (any(!works)) {
        x <- pkg[!works]
        stop(paste0(
            Sys.time(),
            " Package",
            ifelse(length(x) > 1, "s", ""),
            " ",
            paste(x, collapse = ", "),
            " ",
            ifelse(length(x) > 1, "are", "is"),
            " missing. Install ",
            ifelse(length(x) > 1, "them", "it"),
            " with BiocManager::install(",
            ifelse(length(x) > 1, "c(", ""),
            '"',
            paste(x, collapse = '", "'),
            '")',
            ifelse(length(x) > 1, ")", "")
        ))
    }
    return(invisible(NULL))
}
