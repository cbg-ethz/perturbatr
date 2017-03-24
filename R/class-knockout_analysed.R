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


#' Data wrapper for knockout data
#'
#' @rdname knockout_analysed-class
#'
#' @description Class \code{knockout.data} is a wrapper for a
#'   \code{data.table} object
#' containing the knockout data
#'
#' @slot .data the knockout data-set
knockout.analysed <- setClass(
  "knockout.analysed",
  slots     = list(.data="data.table",
                   .inference="character"),
  validity  = function(object)
  {
    if (object@.inference == .inference.types()$MIXED.MODEL)
    {

    }
    else if (object@.inference == .inference.types()$HYPERGEOMETRIC.TEST)
    {

    }
    else if (object@.inference == .inference.types()$T.TEST)
    {

    }
    else stop(paste0("Use one of either ",
                     paste0(.data.types(), collapse="/"),
                     " as data-type."))
    cls <- sort(c("Virus", "Replicate", "Plate",
                  "RowIdx", "ColIdx",
                  "GeneSymbol", "ReadoutType", "Control",
                  "Library", "siRNAIDs", "Screen",
                  "Cell", "ScreenType", "Design",
                  "Entrez", "Readout"))
    if (any(sort(colnames(object@.data)) != cls))
      stop(paste0("Your data needs to have the following colnames:\n",
                  paste0(cls, collapse=", ")))
    return (TRUE)
  }
)
