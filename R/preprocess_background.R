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

#' @noRd
#' @import data.table
#' @importFrom dplyr mutate group_by select
.background.correct <- function(obj, background.column, background.row, method)
{
  ret <-
    dplyr::group_by(obj, Virus, Screen, Library, Replicate, Plate,
                    ScreenType, ReadoutType, ReadoutClass,
                    Design, Cell)

  f <- .summarization.method(method)
  if (!is.null(background.row))
  {
    if (!is.numeric(background.row)) stop("Provide numeric index")
    if (background.row > max(ret$RowIdx))
      stop("Please provide a row index that fits the plate!")
    message(paste("Substracting", method,
                  "of row", background.row, " rom every plate!"))

    ret <- dplyr::mutate(ret, bk = (RowIdx == background.row))
  }
  else if (!is.null(background.column))
  {
    if (!is.numeric(background.column))
      stop("Provide numeric index")
    if (background.column > max(ret$ColIdx))
      stop("Please provide a column index that fits the plate!")
    message(paste("Substracting", method,
                  "of column", background.column, "from every plate!"))

    ret <- dplyr::mutate(ret, bk = (ColIdx == background.column))
  }
  else
  {
    message(paste("Substracting", method, "of all NA genes!"))
    ret <- dplyr::mutate(ret, bk = (is.na(GeneSymbol)))
  }

  ret <-
    dplyr::group_by(ret, Virus, Screen, Library, Replicate, Plate,
                    ScreenType, ReadoutType, ReadoutClass,
                    Design, Cell) %>%
    dplyr::mutate(Readout = .substract.background(Readout, f, bk)) %>%
    ungroup %>%
    dplyr::select(-bk)

  ret
}

#' @noRd
.substract.background <- function(readout, f, background)
{
  ret  <- readout - f(readout[background], na.rm=TRUE)
  ret
}
