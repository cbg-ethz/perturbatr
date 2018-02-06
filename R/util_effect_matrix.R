# perturbR: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2018 Simon Dirmeier
#
# This file is part of perturbR
#
# perturbR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# perturbR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with perturbR. If not, see <http://www.gnu.org/licenses/>.


#' @noRd
#' @param obj  the object to calculate the effect matrices for
#' @param ...  additional parameters
.effect.matrices <- function(obj, ...)
{
  UseMethod(".effect.matrices")
}

#' @noRd
#' @method .effect.matrices perturbation.hm.analysed
#' @import data.table
#' @importFrom dplyr filter select
#' @importFrom tidyr spread
.effect.matrices.perturbation.hm.analysed <- function(obj, ...)
{
  g <- obj@.gene.hits %>%
    dplyr::select(GeneSymbol, Effect) %>%
    .[order(-abs(Effect))]

  pg <- obj@.nested.gene.effects %>%
    dplyr::select(Condition, GeneSymbol, Effect) %>%
    tidyr::spread(Condition, Effect)

  list(gene.effects=g, nested.gene.effects=pg)
}
