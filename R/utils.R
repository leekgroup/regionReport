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
with_wd <- function(dir, expr) {
    wd <- getwd()
    on.exit(setwd(wd))
    setwd(dir)
    eval(expr, envir = parent.frame())
}


#' Attempt to load the namespace of a package and install it if it's missing
#'
#' This function uses requireNamespace to try to load a package. But if it's
#' misssing it will then install it via Bioconductor.
#'
#' @param pkg A single character vector with the name of the package.
#' @param quietly Whether to run requireNamespace and biocLite quietly or not.
#'
load_install <- function(pkg, quietly = TRUE) {
    attemptName <- requireNamespace(pkg, quietly = quietly)
    if(!attemptName) {
        biocLite <- NULL ## To satisfy R CMD check
        
        source('http://bioconductor.org/biocLite.R')
        attemptInstall <- tryCatch(biocLite(pkg, suppressUpdates = quietly),
            warning = function(w) 'failed')
        if(attemptInstall == 'failed') stop(paste('Failed to install', pkg))
        attemptName <- requireNamespace(pkg, quietly = quietly)
    }
    if(attemptName) {
        if(quietly) {
            suppressPackageStartupMessages(library(package = pkg, character.only = TRUE))
        } else {
            library(package = pkg, character.only = TRUE)
        }
    }
    return(invisible(NULL))
}
