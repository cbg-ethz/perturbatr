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


#' @include util_enums.R


#' @aliases show,knockdown.data-method
#' @import data.table
#' @importFrom dplyr select
setMethod(
  "show",
  "knockdown.data",
  function(object)
  {
    cat(paste0("A knockdown data-set\n\n"))
      object@.data[ ,.SD[sample(.N, 2)], by="Virus"] %>%
        dplyr::select(Virus, GeneSymbol, Readout, ScreenType) %>%
        print
  }
)


#' @aliases show,knockdown.lmm.data-method
#' @import data.table
#' @importFrom dplyr select
setMethod(
  "show",
  "knockdown.lmm.data",
  function(object)
  {
    cat(paste0("A knockdown data-set for LMMs\n\n"))
    object@.data[ ,.SD[sample(.N, 2)], by="Virus"] %>%
      dplyr::select(Virus, GeneSymbol, Readout, ScreenType, Weight) %>%
      print
  }
)


#' @aliases show,knockdown.lmm.analysed-method
#' @import data.table
#' @importFrom dplyr select left_join
#' @importFrom tidyr spread
setMethod(
  "show",
  "knockdown.lmm.analysed",
  function(object)
  {
    cat(paste0("An LMM-analyed knockdown data-set\n\n"))
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
  "knockdown.hyper.analysed",
  function(object)
  {
    cat(paste0("An analyed knockdown data-set using an ",
               "iterative hypergeometric test\n\n"))
    data.table::setorder(object@.gene.hits, -HitRatio, MinQval)
    object@.gene.hits[ ,head(.SD, 2L), by="Virus"] %>%
      dplyr::select(Virus, GeneSymbol, MeanEffect, HitRatio, MinQval) %>%
      print
  }
)

