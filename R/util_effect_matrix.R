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
effect.matrices <- function(obj)
{
  UseMethod("effect.matrices")
}

#' @noRd
#' @method .effect.matrices HMAnalysedPerturbationData
#' @import tibble
#' @importFrom dplyr filter select
#' @importFrom tidyr spread
effect.matrices.HMAnalysedPerturbationData <- function(obj)
{
  g <- geneHits(obj) %>%
    dplyr::select(GeneSymbol, Effect) %>%
    dplyr::arrange(desc(abs(Effect)))

  pg <- nestedGeneEffects(obj) %>%
    dplyr::select(Condition, GeneSymbol, Effect) %>%
    tidyr::spread(Condition, Effect)

  list(gene.effects=g, nested.gene.effects=pg)
}
