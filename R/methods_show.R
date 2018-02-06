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


#' @aliases show,perturbation.data-method
#' @import data.table
#' @importFrom dplyr select
setMethod(
	"show",
	"perturbation.raw.data",
	function(object)
	{
		cat(paste0("A raw perturbation data-set\n\n"))
		object@.data[ ,.SD[sample(.N, 2)], by="Condition"] %>%
			dplyr::select(Condition, GeneSymbol, Readout) %>%
			print
	}
)

#' @aliases show,perturbation.data-method
#' @import data.table
#' @importFrom dplyr select
setMethod(
  "show",
  "perturbation.data",
  function(object)
  {
    cat(paste0("A perturbation data-set\n\n"))
      object@.data[ ,.SD[sample(.N, 2)], by="Condition"] %>%
        dplyr::select(Condition, GeneSymbol, Readout) %>%
        print
  }
)


#' @aliases show,perturbation.hm.data-method
#' @import data.table
#' @importFrom dplyr select
setMethod(
  "show",
  "perturbation.hm.data",
  function(object)
  {
    cat(paste0("A perturbation data-set for HM\n\n"))
    object@.data[ ,.SD[sample(.N, 2)], by="Condition"] %>%
      dplyr::select(Condition, GeneSymbol, Readout, Weight) %>%
      print
  }
)


#' @aliases show,perturbation.hm.analysed-method
#' @import data.table
#' @importFrom dplyr select left_join
#' @importFrom tidyr spread
setMethod(
  "show",
  "perturbation.hm.analysed",
  function(object)
  {
    cat(paste0("A perturbation data-set analysed using a hierachical model\n\n"))
    gps <- object@.nested.gene.effects %>%
      dplyr::select(GeneSymbol, Condition, Effect) %>%
      tidyr::spread(Condition, Effect)
    ges <- object@.gene.effects %>%
      dplyr::select(GeneSymbol, Effect, Qval)
    mer <- dplyr::left_join(ges, gps, by="GeneSymbol")
    print(data.table::as.data.table(mer))
  }
)

#' @aliases show,perturbation.analysed-method
#' @import data.table
#' @importFrom dplyr select
setMethod(
  "show",
  "perturbation.analysed",
  function(object)
  {
    cat(paste0("A perturbation data-set analysed using a ", object@.inference, "\n\n"))
    print(object@.gene.hits)
  }
)
