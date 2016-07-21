#' @noRd
#' @export
#' @import data.table
print.svd.plates    <- function(x, ...) .print.svd.list(x, ...)

#' @noRd
#' @export
#' @import data.table
print.svd.replicates <- function(x, ...) .print.svd.list(x, ...)

#' @noRd
#' @import data.table
.print.svd.list <-
  function
(
  x,
  ...
)
{
  if (length(x) != 0)
  {
    cat(paste("List of", length(x), "replicates! Showing structure of plates.\n\n"))
    return(str(x[[1]]))
  }
  else return("Empty list!\n")
}

#' @noRd
#' @export
#' @import data.table
#' @importFrom dplyr group_by summarize as.tbl
print.svd.data <-
function
(
  x,
  ...
)
{
  if (nrow(x) != 0)
  {
    ret <-
      dplyr::group_by(x, Virus, Screen, Replicate,
                      InfectionType, ReadoutType, Design, Cell, Library) %>%
      dplyr::summarize(Plates=max(Plate), Wells=(max(ColIdx)*max(RowIdx))) %>%
      ungroup
    cat("Printing data overview!\n")
    print(data.table:::print.data.table(ret))
  }
  cat("\n\nPrinting data!\n")
  print(data.table:::print.data.table(x))
}

#' @noRd
#' @export
#' @import data.table
#' @importFrom dplyr group_by summarize filter
print.svd.raw <-
function
(
  x,
  ...
)
{
  cat("Filtering by 'readout', hiding 'viability'!\n\n")
  ret <-
    dplyr::filter(x, ReadoutClass=="Readout")
  return(print.svd.data(ret))
}

#' @noRd
#' @export
#' @import data.table
#' @importFrom dplyr group_by summarize as.tbl
print.svd.analysed.hyper <-
function
(
  x,
  ...
)
{
  ret <-
    dplyr::group_by(x, Virus, Screen, InfectionType, ReadoutType,
                    Design, Cell, Library, GeneSymbol) %>%
    dplyr::summarize(siRNAIDs=n()) %>%
    ungroup
  cat("Printing data overview!\n")
  print(data.table:::print.data.table(ret))
  cat("\n\nPrinting data!\n")
  print(data.table:::print.data.table(x))
}

#' @noRd
#' @export
#' @import data.table
print.svd.prioritized.hyper <- function(x, ...) .print.svd.prioritized(x, ...)

#' @noRd
#' @export
#' @import data.table
print.svd.prioritized.tt <- function(x, ...) .print.svd.prioritized(x, ...)

#' @noRd
#' @import data.table
#' @importFrom dplyr group_by filter filter
.print.svd.prioritized <-
function
(
  x,
  ...
)
{
  ret <-
    x[, .SD[order(abs(MeanEffect), decreasing=T)[1:25]],
      by=c("Virus", "Screen", "ReadoutType", "InfectionType",
           "Library", "Design", "Cell")] %>%
    dplyr::filter(!is.na(GeneSymbol))
  cat("Printing data overview!\n")
  print(data.table:::print.data.table(ret))
  return
}

#' @noRd
#' @export
#' @import data.table
print.svd.quality <-
function
(
  x,
  ...
)
{
  q <- x$quality
  cat("Plate quality:\n")
  print(data.table:::print.data.table(q$plate.quality))
  cat("Replicate quality:\n")
  print(data.table:::print.data.table(q$replicate.quality))
  cat("Screen quality:\n")
  print(data.table:::print.data.table(q$screen.quality))
}
