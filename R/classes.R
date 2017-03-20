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
#' @rdname knockout_data-class
#'
#' @description Class \code{knockout.data} is a wrapper for a \code{data.table} object
#' containing the knockout data
#'
#' @slot .data the knockout data-set
knockout.data <- setClass(
  "knockout.data",
  slots     = list(.data="data.table"),
  validity  = function(object)
  {
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

#' Data wrapper for unnormalized knockout data
#'
#' @rdname knockout_raw_data-class
#'
#' @description Class \code{knockout.data.raw} is a wrapper for a \code{data.table} object
#' containing the unnormalized knockout data
#'
#' @slot .data the unnormalized knockout data-set
knockout.data.raw <- setClass(
  "knockout.data.raw",
  slots     = list(.data="data.table"),
  validity  = function(object)
  {
    cls <- sort(c("Virus", "Replicate", "Plate",
                  "RowIdx", "ColIdx",
                  "GeneSymbol", "ReadoutClass", "ReadoutType", "Control",
                  "Library", "siRNAIDs", "Screen",
                  "Cell", "ScreenType", "Design",
                  "Entrez", "Readout", "NumCells"))
    if (any(sort(colnames(object@.data)) != cls))
      stop(paste0("Your data needs to have the following colnames:\n",
                  paste0(cls, collapse=", ")))
    return (TRUE)
  }
)

