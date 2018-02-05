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
# along with perturbR If not, see <http://www.gnu.org/licenses/>.


#' @include util_enums.R


#' @noRd
#' @slot .data the perturbation data-set
setClass(
  "perturbation.data",
  contains = "VIRTUAL",
  slots    = list(.data="data.table"),
)


#' @title Data wrapper for raw perturbation data
#'
#' @description Class \code{perturbation.raw.data} is a wrapper for a
#' \code{data.table} object containing the raw data-set
#'
#' @name Raw-data
#' @rdname perturbation_raw_data-class
#'
setClass(
  "perturbation.raw.data",
  contains  = "perturbation.data",
  validity  = function(object)
  {
    cls <- c("Virus", "Replicate", "Plate", "RowIdx", "ColIdx",
             "GeneSymbol", "ReadoutType", "Control", "Library",
             "siRNAIDs", "Screen", "Cell", "ScreenType", "Design",
             "Entrez", "Readout", "ReadoutClass", "NumCells")
    if (any(sort(colnames(object@.data)) != sort(cls)))
      stop(paste0("Your data needs to have the following colnames:\n",
                  paste0(cls, collapse=", ")))

  return (TRUE)
  }
)

#' @title Data wrapper for normalized perturbation data
#'
#' @description Class \code{perturbation.normalized.data} is a wrapper for a
#' \code{data.table} object containing the normalized data-set
#'
#' @name Normalized-data
#' @rdname perturbation_normalized_data-class
#'
setClass(
  "perturbation.normalized.data",
  contains  = "perturbation.data",
  validity  = function(object)
  {
      cls <- c("Virus", "Replicate", "Plate", "RowIdx", "ColIdx",
               "GeneSymbol", "ReadoutType", "Control", "Library",
               "siRNAIDs", "Screen", "Cell", "ScreenType", "Design",
               "Entrez", "Readout")
      if (any(sort(colnames(object@.data)) != sort(cls)))
        stop(paste0("Your data needs to have the following colnames:\n",
                    paste0(cls, collapse=", ")))
    return (TRUE)
  }
)

#' @title Data wrapper for linear-mixed-model perturbation data
#'
#' @description Class \code{perturbation.lmm.data} is a wrapper for a
#' \code{data.table} object containing the normalized data-set
#'
#' @name LMM-data
#' @rdname perturbation_lmm_data-class
#'
setClass(
  "perturbation.lmm.data",
  contains  = "perturbation.data",
  validity  = function(object)
  {
    cls <- c("Virus", "GeneSymbol", "ReadoutType", "Control", "Weight",
             "Cell", "ScreenType", "Design", "Entrez", "Readout", "VG")
    if (any(!(colnames(object@.data) %in% sort(cls))))
      stop(paste0("Your data needs to have the following colnames:\n",
                  paste0(cls, collapse=", ")))
    return (TRUE)
  }
)
