#' Extend the results from a <code>svd.prioritized</code> object by network diffusion.
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
#' @param path  path to the network
#' @param ...  additional parameters
diffuse <- function(obj, method=c("neighbors", "mrw"), path, ...) UseMethod("diffuse")

#' @noRd
#' @export
#' @import data.table
diffuse.svd.prioritized.pmm <- function(obj, method=c("neighbors", "mrw"), path, ...)
{
  if (!file.exists(path)) stop(paste("Can't find: ", path, "! Yieks!", sep=""))
  graph <- .read.graph(file)
  res  <- .diffuse.lmm(obj, match.arg(method), graph)
  class(res) <- c("svd.diffused.pmm","svd.diffused", class(res))
  invisible(res)
}

#' @noRd
#' @import data.table
.diffuse.lmm <- function(obj, method, graph)
{

}


