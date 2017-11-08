# knockdown: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2015 - 2016 Simon Dirmeier
#
# This file is part of knockdown
#
# knockdown is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# knockdown is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with knockdown. If not, see <http://www.gnu.org/licenses/>.


#' Get the plate matrix (plus control indexes) from an object
#'
#' @noRd
#' @param obj  the object for which you want to have the readout plate.matrix
#' @param ...  additional params
.plate.matrix <- function(obj, ...) UseMethod(".plate.matrix")

#' @noRd
#' @export
#' @method .plate.matrix knockdown.plate
#' @import data.table
#' @importFrom dplyr filter
#' @importFrom dplyr select
.plate.matrix.knockdown.plate <- function(obj, ...)
{

  # highly inefficient method!!
  # :(

  col.size <- max(obj@.data$ColIdx)
  row.size <- max(obj@.data$RowIdx)
  m <- idx <- genes <- matrix(0, row.size, col.size)
  # highly inefficient: Rcpp?
  for (i in 1:row.size)
  {
    row <- dplyr::filter(obj@.data, RowIdx==i)
    for (j in 1:col.size)
    {
      col <- dplyr::filter(row, ColIdx==j)
      if (nrow(col) == 0)
      {
        m[i, j]    <- NA_real_
        idx[i, j]  <- 0
        genes[i,j] <- NA_character_
      }
      else
      {
        m[i, j]     <- dplyr::select(col, Readout) %>% unlist
        idx[i, j]   <- dplyr::select(col, Control)  %>% unlist
        genes[i, j] <- dplyr::select(col, GeneSymbol) %>% unlist
      }
    }
  }
  res <- list(readout=m, idx=idx, genes=genes)
  res
}
