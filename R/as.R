#' Converts an object to an svd.data object
#'
#' @export
#' @import data.table
#' @param obj  the object to be converted
as.svd.data <- function(obj, ...) UseMethod("as.svd.data")

#' @export
as.svd.data.data.table <-
function
(
  obj,
  ...
)
{
  class(obj) <- c("svd.data", class(obj))
  invisible(obj)
}
