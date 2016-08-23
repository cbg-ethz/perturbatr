#' Fit an LMM to the data and calculate local false discovery rates.
#'
#' @export
#'
#' @import data.table
#'
#' @param obj  an svd.data object
#' @param method  method that should be used for diffusion
#' \itemize{
#'   \item{neighbors }{ just looks at the neighbors :)}
#'   \item{neighbors }{ just looks at the neighbors :)}
#' }
#' @param ...  additional parameters
diffuse <- function(obj, method=c("neighbors", "mrw"), ...) UseMethod("diffuse")

#' @noRd
#' @export
#' @import data.table
diffuse.svd.prioritized.pmm <- function(obj, method=c("neighbors", "mrw"), ...)
{
  res  <- .diffuse.lmm(obj, match.arg(method))
  class(res) <- c("svd.analysed.pmm","svd.analysed", class(res))
  invisible(res)
}




