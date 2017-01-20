
#' Plot a table to the grid
#'
#' @export
#' @param obj  the object you want
show.table <- function(obj, nrow=10, colored=2, main="", ...)
{
  stop{"todo"}
  UseMethod("show.table")
}

#' @export
#' @import xtable
#' @import
#' @method show.table data.table
show.table.data.table <- function(obj, nrow=10, colored=2, main="", ...)
{
  print(xtable(t, align=c("lllll"), digits=5, label="tab:lmm-hits",
               caption="Benjamini-Hochberg corrected P-values of 10 strongest hits."),
        booktabs=T, comment=FALSE, sanitize.colnames.function=bold)
}
