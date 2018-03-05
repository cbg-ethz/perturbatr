# perturbatr: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2018 Simon Dirmeier
#
# This file is part of perturbatr
#
# perturbatr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# perturbatr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with perturbatr. If not, see <http://www.gnu.org/licenses/>.


#' @noRd
#' @import data.table
#' @importFrom dplyr mutate group_by
.remove.outliers <- function(obj, rm.outlier.wells)
{
  stopifnot(length(rm.outlier.wells) == 2)
  if (any(rm.outlier.wells > 1)) stop("Some outlier quantiles >1")
  if (any(rm.outlier.wells < 0)) stop("Some outlier quantiles <0")
  .rm.wells.quantile(obj, outlier.well.range),
}

#' @noRd
.rm.wells.quantile <- function(obj, probs)
{
  stopifnot(all(probs <= 1), all(probs >= 0))
  message(paste("Removing wells based on ", probs[1],
                "% and ", probs[2], "% quantile of cell number!"))
  obj <- dplyr::group_by(obj, Condition, Screen, Library,
                         ScreenType, ReadoutType, ReadoutClass,
                         Cell, Design) %>%
    dplyr::mutate(Readout=.rmwqs(Readout, NumCells, probs)) %>%
    ungroup
  obj
}

#' @noRd
#' @importFrom stats quantile
.rmwqs <- function(read, num, probs)
{
  re <- read
  if (!all(is.na(num)))
  {
    quant <- unname(stats::quantile(num, na.rm=TRUE, probs=probs))
    message(paste0("\t..removing x<", quant[1], " | ", quant[2], "<x wells!"))
    re[num < quant[1] | quant[2] < num] <- NA_real_
  }
  re
}
