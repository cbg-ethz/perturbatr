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

#' Calculcate the PMM-prioritized effect matrices for gene and pathogen-gene matrices
#'
#' @export
#' @import data.table
#' @param obj  the object to calculate the effect matrices for
#' @param ...  additional parameters
effect.matrices <- function(obj, ...)
{
  UseMethod("effect.matrices")
}

#' @noRd
#' @export
#' @import data.table
#' @importFrom dplyr filter select
effect.matrices.svd.prioritized.pmm <- function(obj, ...)
{

  gene.top <- obj$gene.hits %>%
    dplyr::select(GeneSymbol, Effect) %>%
    .[order(-abs(Effect))] %>%
    .[, .SD[1:min(25,.N)]]
  pgs <- obj$fit$gene.pathogen.effects %>%
    dplyr::filter(GeneSymbol %in% gene.top$GeneSymbol) %>%
    dplyr::select(Virus, GeneSymbol, Effect) %>%
    tidyr::spread(Virus, Effect)
  res <- list(gene.effects=gene.top, gene.pathogen.hits=pgs)
  class(res) <- "svd.prioritized.pmm.effect.matrices"
  res
}
