#' Deprecated functions in package \sQuote{regionReport}
#'
#' These functions are provided for compatibility with older versions
#' of \sQuote{regionReport} only, and will be defunct at the next release.
#'
#' @details 
#'
#' The following functions are deprecated and will be made defunct; use
#'  the replacement indicated below:
#'  \itemize{
#'
#'    \item{plot_cluster: \code{\link{plotCluster}}}
#'    \item{render_report: \code{\link{render_report}}}
#'
#'  }
#'
#' @name regionReport-deprecated
#' @export
derfinder_report <- function() { .Deprecated('derfinderReport')}

#' @rdname regionReport-deprecated
#' @export
render_report <- function() { .Deprecated('renderReport')}
