#' @noRd
#' @export
#' @import data.table
print.svd.plates    <- function(x, ...) .print.svd.list(x, ...)

#' @noRd
#' @export
#' @import data.table
print.svd.replicates <- function(x, ...) .print.svd.list(x, ...)

#' @noRd
#' @export
#' @import data.table
#' @importFrom utils str
.print.svd.list <-
  function
(
  x,
  ...
)
{
  if (base::length(x) != 0)
  {
    cat(base::paste("List of", base::length(x),
                    "replicates! Showing structure of plates.\n\n"))
    return(utils::str(x[[1]]))
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
    base::cat("Printing data overview!\n")
    base::print(data.table:::print.data.table(ret))
  }
  base::cat("\n\nPrinting data!\n")
  base::print(data.table:::print.data.table(x))
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
  base::cat("Printing data overview!\n")
  base::print(data.table:::print.data.table(ret))
  cat("\n\nPrinting data!\n")
  base::print(data.table:::print.data.table(x))
}

#' @noRd
#' @export
#' @import data.table
#' @importFrom dplyr group_by summarize as.tbl
print.svd.analysed.pmm <-
function
(
  x,
  ...
)
{
  ret <- x$gene.pathogen.effects %>%
    dplyr::group_by(Virus) %>%
    dplyr::summarize(GenesCount=n()) %>%
    ungroup
  cat("Printing data overview for gene-pathogen effects!\n")
  print(data.table:::print.data.table(ret))
  cat("\n\nPrinting data for gene-pathogen effects!\n")
  print(data.table:::print.data.table(x$gene.pathogen.effects))
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
#' @export
#' @import data.table
#' @importFrom dplyr filter
print.svd.prioritized.pmm <-
function
(
  x,
 ...
)
{
  ret <-
    x$gene.pathogen.effect.hits[, .SD[order(abs(Effect), decreasing=T)[1:25]],
      by=c("Virus")] %>%
    dplyr::filter(!is.na(GeneSymbol))
  cat("Printing data overview for virus-gene effects!\n")
  print(data.table:::print.data.table(ret))
  ret <-
    x$gene.effect.hits[, .SD[order(abs(Effect), decreasing=T)[1:25]]] %>%
    dplyr::filter(!is.na(GeneSymbol))
  cat("Printing data overview for gene effects!\n")
  print(data.table:::print.data.table(ret))
  return
}

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
  cat("Printing data overview for!\n")
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
