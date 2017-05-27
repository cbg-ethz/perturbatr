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


#' Get the readout matrix (plus control indexes) from an object
#'
#' @export
#'
#' @param obj  the object for which you want to have the readout matrix
#' @param ...  additional params
readout.matrix <- function(obj, ...) UseMethod("readout.matrix")

#' @export
#' @method readout.matrix knockout.plate
#' @import data.table
#' @importFrom dplyr filter
#' @importFrom dplyr select
readout.matrix.knockout.plate <- function(obj, ...)
{
  col.size <- max(obj@.data$ColIdx)
  row.size <- max(obj@.data$RowIdx)
  m <- idx <- genes <- matrix(0, row.size, col.size)
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
