#' @noRd
#' @import data.table
#' @importFrom dplyr mutate
.log.norm <- function(obj)
{
  message("Calculating log!")
  obj <- dplyr::mutate(obj, Readout = log(Readout  + 0.00001))
  invisible(obj)
}
