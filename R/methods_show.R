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


#' @include util_enums.R


#' @aliases show,PerturbationData-method
#' @import tibble
#' @importFrom dplyr select group_by sample_n
setMethod(
  "show",
  "PerturbationData",
  function(object)
  {
    cat("A perturbation data set\n\n")
    dat <- dplyr::group_by(dataSet(object), Condition) %>%
      dplyr::sample_n(2, replace=TRUE) %>%
      dplyr::select(Condition, GeneSymbol, Readout)
    print(tibble::as.tibble(dat))
  }
)


#' @aliases show,HMAnalysedPerturbationData-method
#' @import tibble
#' @importFrom dplyr select left_join
#' @importFrom tidyr spread
setMethod(
  "show",
  "HMAnalysedPerturbationData",
  function(object)
  {
    cat(paste0(
      "A perturbation data-set analysed using a hierachical model\n\n"))
    gps <- nestedGeneEffects(object) %>%
      dplyr::select(GeneSymbol, Condition, Effect) %>%
      tidyr::spread(Condition, Effect)
    ges <- geneEffects(object) %>%
      dplyr::select(GeneSymbol, Effect, Qval)
    mer <- dplyr::left_join(ges, gps, by="GeneSymbol")
    print(tibble::as.tibble(mer))
  }
)


#' @aliases show,NetworkAnalysedPerturbationData-method
#' @import tibble
#' @importFrom dplyr select arrange desc
setMethod(
  "show",
  "NetworkAnalysedPerturbationData",
  function(object)
  {
    cat(paste0(
      "A perturbation data-set analysed usingnetwork diffuson\n\n"))
    gps <- geneEffects(object) %>%
      dplyr::select(GeneSymbol, Effect, DiffusionEffect)
    gps <- gps[order(-gps$DiffusionEffect),]
    print(tibble::as.tibble(gps))
  }
)
