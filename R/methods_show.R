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


#' @include util_enums.R


#' @aliases show,perturbation.data-method
#' @import data.table
#' @importFrom dplyr select
setMethod(
  "show",
  "perturbation.data",
  function(object)
  {
    cat(paste0("A perturbation data-set\n\n"))
      object@.data[ ,.SD[sample(.N, 2)], by="Virus"] %>%
        dplyr::select(Virus, GeneSymbol, Readout, ScreenType) %>%
        print
  }
)


#' @aliases show,perturbation.lmm.data-method
#' @import data.table
#' @importFrom dplyr select
setMethod(
  "show",
  "perturbation.lmm.data",
  function(object)
  {
    cat(paste0("A perturbation data-set for LMMs\n\n"))
    object@.data[ ,.SD[sample(.N, 2)], by="Virus"] %>%
      dplyr::select(Virus, GeneSymbol, Readout, ScreenType, Weight) %>%
      print
  }
)


#' @aliases show,perturbation.lmm.analysed-method
#' @import data.table
#' @importFrom dplyr select left_join
#' @importFrom tidyr spread
setMethod(
  "show",
  "perturbation.lmm.analysed",
  function(object)
  {
    cat(paste0("An LMM-analyed perturbation data-set\n\n"))
    gps <- object@.gene.pathogen.effects %>%
      dplyr::select(GeneSymbol, Virus, Effect) %>%
      tidyr::spread(Virus, Effect)
    ges <- object@.gene.effects %>%
      dplyr::select(GeneSymbol, Effect, Qval)
    mer <- dplyr::left_join(ges, gps, by="GeneSymbol")
    print(data.table::as.data.table(mer))
  }
)


#' @import data.table
setMethod(
  "show",
  "perturbation.hyper.analysed",
  function(object)
  {
    cat(paste0("An analyed perturbation data-set using an ",
               "iterative hypergeometric test\n\n"))
    data.table::setorder(object@.gene.hits, -HitRatio, MinQval)
    object@.gene.hits[ ,head(.SD, 2L), by="Virus"] %>%
      dplyr::select(Virus, GeneSymbol, MeanEffect, HitRatio, MinQval) %>%
      print
  }
)

