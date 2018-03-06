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
# along with perturbatr If not, see <http://www.gnu.org/licenses/>.


#' @include util_enums.R
#' @include methods_getters.R


#' @noRd
.requiredDataCols  <- function()
{
    c("Condition", "Replicate", "GeneSymbol", "Perturbation",
       "Readout", "Control")
}


#' @title Data wrapper for a data set of a perturbation screen.
#'
#' @description Class \code{PerturbationData} wraps a data set derived
#'  from a genetic perturbation screen, e.g. using RNA interference or CRISPR.
#'  Class \code{PerturbationData} exposes getters for its members of the same
#'  name, but no setters, because the data should be treated as
#'  constant once set. The easiest way to construct a \code{PerturbationData}
#'  object is by first creating a \code{data.frame} and then calling \code{as}.
#'  See the examples to construct an object.
#'
#' @name PerturbationData-class
#' @rdname PerturbationData-class
#' @exportClass PerturbationData
#'
#' @slot dataSet  the data set as a \code{data.table} object
#' @examples
#'   df <- data.frame(Condition    = c("V1", "V2", "V3"),
#'                    Replicate    = c(1, 1, 1),
#'                    GeneSymbol   = c("TP52", "NegCtrl", "PosCtrl"),
#'                    Perturbation = c("P1", "P2", "P3"),
#'                    Readout      = c(123, 121, 12),
#'                    Control      = c(0, -1, 1))
#'   methods::as(df, "PerturbationData")
#'
setClass(
    "PerturbationData",
    slots    = list(dataSet="data.table"),
    validity = function(object)
    {
        check.columns(object@dataSet, .requiredDataCols())
        return(TRUE)
    }
)


#' @rdname dataSet-methods
#' @aliases dataSet,PerturbationData-method
#' @import data.table
setMethod(
    "dataSet",
    signature = signature(obj="PerturbationData"),
    function(obj) obj@dataSet)
