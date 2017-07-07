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


#' @noRd
#' @import data.table
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
.loess <- function(obj)
{
  message("Calculating LOESS model!")
  loessed.data <-
    dplyr::group_by(obj, Virus, Screen, Library,
                    ReadoutType, ScreenType, ReadoutClass,
                    Design, Cell,
                    Replicate, Plate) %>%
    dplyr::mutate(Readout = .loess.plate(NumCells, Readout)) %>%
    ungroup
  invisible(loessed.data)
}

#' @noRd
#' @importFrom stats lowess
.loess.plate <- function(n.cells, readout)
{
  good.idxs <- (!is.na(n.cells) & !is.na(readout))
  sorted.n.cells     <- sort.int(n.cells[good.idxs], index.return = TRUE)
  sorted.n.cells.idx <- sort.int(sorted.n.cells$ix,  index.return = TRUE)
  loessed.readout <- rep(NA_real_, length(readout))
  if (all(is.na(readout)) || all(is.na(n.cells)))
  {
    loessed.readout <- readout
  }
  else
  {
    tryCatch({
      res <- stats::lowess(n.cells[good.idxs], readout[good.idxs])
      loessed.readout[good.idxs] <-
        readout[good.idxs] - res$y[sorted.n.cells.idx$ix] },
      warning = function(war) { message(paste("WARNING: ", war, "\n")); },
      error = function(err)   { message(paste("ERROR: ", err,  "\n")); }
    )
  }
  if (all(is.na(loessed.readout))) loessed.readout <- readout
  invisible(loessed.readout)
}
