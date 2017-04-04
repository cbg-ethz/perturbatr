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


#' @include util_enums.R


#' @aliases show,knockout.data-method
#' @import data.table
#' @importFrom dplyr select
setMethod(
  "show",
  "knockout.data",
  function(object)
  {
    cat(paste0("A knockout data-set\n\n"))
      object@.data[ ,.SD[sample(.N, 2)], by="Virus"] %>%
        dplyr::select(Virus, GeneSymbol, Readout, Library,
                      ReadoutType, Screen, Cell, ScreenType, Design) %>%
        print
  }
)


#' @aliases show,knockout.lmm.data-method
#' @import data.table
#' @importFrom dplyr select
setMethod(
  "show",
  "knockout.lmm.data",
  function(object)
  {
    cat(paste0("A knockout data-set for LMMs\n\n"))
    object@.data[ ,.SD[sample(.N, 2)], by="Virus"] %>%
      dplyr::select(Virus, GeneSymbol, Readout, Weight,
                    ReadoutType, Cell, ScreenType, Design) %>%
      print
  }
)


#' @aliases show,knockout.lmm.analysed-method
#' @import data.table
#' @importFrom dplyr select left_join
#' @importFrom tidyr spread
setMethod(
  "show",
  "knockout.lmm.analysed",
  function(object)
  {
    cat(paste0("An LMM-analyed knockout data-set\n\n"))
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
  "knockout.hyper.analysed",
  function(object)
  {
    cat(paste0("An analyed knockout data-set using an ",
               "iterative hypergeometric test\n\n"))
    data.table::setorder(object@.gene.hits, -HitRatio, MinQval)
    object@.gene.hits[ ,head(.SD, 2L), by="Virus"] %>%
      dplyr::select(Virus, GeneSymbol, MeanEffect, HitRatio, MinQval) %>%
      print
  }
)

