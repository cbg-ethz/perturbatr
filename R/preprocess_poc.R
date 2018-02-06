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
.poc <- function(obj, method, ctrl.gene)
{
  ctrl.gene <- tolower(ctrl.gene)
  if (!is.na(ctrl.gene))
    message(paste0("Calculating POC on gene ", ctrl.gene, "!"))
  else message("Calculating POC on negative controls!")
  f <- .summarization.method(method)
  ret <-
    dplyr::group_by(obj, Condition, Screen, Library,
                    ReadoutType, ScreenType, ReadoutClass,
                    Design, Cell, Replicate, Plate) %>%
    dplyr::mutate(Readout=.poc.grp(Readout, ctrl.gene, f, Control, GeneSymbol))
  ret
}

#' @noRd
.poc.grp <- function(read, ctrl.gene, f, ctrl, genes)
{
  idxs <- which(ctrl == -1)
  if (!is.na(ctrl.gene)) idxs <- which(ctrl == -1 & ctrl.gene == genes)
  if (length(idxs) == 0) stop("No normalization genes found!")
  ret <- (read / f(read[idxs])) * 100
  ret
}
