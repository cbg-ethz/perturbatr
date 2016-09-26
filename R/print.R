# knockout: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2015 - 2016 Simon Dirmeier
#
# This file is part of knockout
#
# knockout is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# knockout is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with knockout. If not, see <http://www.gnu.org/licenses/>.

#' @export
#' @import data.table
#' @method print svd.plates
print.svd.plates <- function(x, ...) .print.svd.list(x, ...)

#' @export
#' @import data.table
#' @method print svd.replicates
print.svd.replicates <- function(x, ...) .print.svd.list(x, ...)

#' @noRd
#' @import data.table
#' @importFrom utils str
.print.svd.list <- function(x, ...)
{
  if (base::length(x) != 0)
  {
    cat(base::paste("List of", base::length(x),
                    "replicates! Showing structure of plates.\n\n"))
    return(utils::str(x[[1]]))
  }
  else return("Empty list!\n")
}

#' @export
#' @import data.table
#' @importFrom dplyr group_by summarize as.tbl
#' @method print svd.data
print.svd.data <- function(x, ...)
{
  if (nrow(x) != 0)
  {
    ret <-
      dplyr::group_by(x, Virus, Screen, Replicate,
                      InfectionType, ReadoutType, Design, Cell, Library) %>%
      dplyr::summarize(Plates=max(Plate), Wells=(max(ColIdx)*max(RowIdx))) %>%
      ungroup
    base::cat("Printing data overview!\n")
    base::print(data.table::as.data.table(ret))
  }
  base::cat("\n\nPrinting data!\n")
  base::print(data.table::as.data.table(x))
}

#' @export
#' @import data.table
#' @importFrom dplyr group_by summarize filter
#' @method print svd.raw
print.svd.raw <- function(x, ...)
{
  cat("Filtering by 'readout', hiding 'viability'!\n\n")
  ret <-
    dplyr::filter(x, ReadoutClass=="Readout")
  return(print.svd.data(ret))
}

#' @export
#' @import data.table
#' @importFrom dplyr group_by summarize as.tbl
#' @method print svd.analysed.hyper
print.svd.analysed.hyper <- function(x,  ...)
{
  ret <-
    dplyr::group_by(x, Virus, Screen, InfectionType, ReadoutType,
                    Design, Cell, Library, GeneSymbol) %>%
    dplyr::summarize(siRNAIDs=n()) %>%
    ungroup
  base::cat("Printing data overview!\n")
  base::print(data.table::as.data.table(ret))
  cat("\n\nPrinting data!\n")
  base::print(data.table::as.data.table(x))
}



#' @export
#' @import data.table
#' @importFrom dplyr group_by summarize as.tbl
#' @method print svd.analysed.pmm
print.svd.analysed.pmm <- function(x, ...)
{
  ret <- x$gene.pathogen.effects %>%
    dplyr::group_by(Virus) %>%
    dplyr::summarize(GenesCount=n()) %>%
    ungroup
  cat("Printing data overview for gene-pathogen effects!\n")
  print(data.table::as.data.table(ret))
  cat("\n\nPrinting data for gene-pathogen effects!\n")
  print(data.table::as.data.table(x$gene.pathogen.effects))
}

#' @export
#' @import data.table
#' @method print svd.prioritized.hyper
print.svd.prioritized.hyper <- function(x, ...) .print.svd.prioritized(x, ...)


#' @export
#' @import data.table
#' @method print svd.prioritized.tt
print.svd.prioritized.tt <- function(x, ...) .print.svd.prioritized(x, ...)

#' @export
#' @import data.table
#' @importFrom dplyr filter
#' @method print svd.prioritized.pmm
print.svd.prioritized.pmm <- function(x, ...)
{
  ret <-
    x$gene.pathogen.effect.hits[, .SD[order(abs(Effect), decreasing=T)[1:25]],
      by=c("Virus")] %>%
    dplyr::filter(!is.na(GeneSymbol))
  cat("Printing data overview for virus-gene effects!\n")
  print(data.table::as.data.table(ret))
  ret <-
    x$gene.effect.hits[, .SD[order(abs(Effect), decreasing=T)[1:25]]] %>%
    dplyr::filter(!is.na(GeneSymbol))
  cat("Printing data overview for gene effects!\n")
  print(data.table::as.data.table(ret))
  return
}

#' @noRd
#' @import data.table
#' @importFrom dplyr group_by filter filter
.print.svd.prioritized <- function(x, ...)
{
  ret <-
    x[, .SD[order(abs(MeanEffect), decreasing=T)[1:25]],
      by=c("Virus", "Screen", "ReadoutType", "InfectionType",
           "Library", "Design", "Cell")] %>%
    dplyr::filter(!is.na(GeneSymbol))
  cat("Printing data overview for!\n")
  print(data.table::as.data.table(ret))
  return
}

#' @export
#' @import data.table
#' @method print svd.quality
print.svd.quality <- function(x, ...)
{
  q <- x$quality
  cat("Plate quality:\n")
  print(data.table::as.data.table(q$plate.quality))
  cat("Replicate quality:\n")
  print(data.table::as.data.table(q$replicate.quality))
  cat("Screen quality:\n")
  print(data.table::as.data.table(q$screen.quality))
}
