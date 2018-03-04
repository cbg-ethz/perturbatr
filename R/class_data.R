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


#' @noRd
.requiredDataCols  <- function()
{
    c("Condition", "Replicate", "GeneSymbol", "Perturbation", "Readout")
}


#' @noRd
.requiredHMCols   <- function()
{
    c("Condition", "GeneSymbol", "Weight", "Readout")
}


#' @title Data wrapper for a data set of a perturbation screen.
#'
#' @description Class \code{PerturbationData} wraps a data set derived
#'  from a genetic perturbation screen, e.g. using RNA interference or CRISPR.
#'
#' @name PerturbationData-class
#' @rdname PerturbationData-class
#' @exportClass PerturbationData
#'
#' @slot dataSet  the data set as a \code{data.table} object
#' @slot dataType  the type that describes your data best, for instance
#'  \code{raw} or \code{normalized}
#'
setClass(
    "PerturbationData",
    slots    = list(dataSet="data.table", dataType="character")
    validity = function(object)
    {
        if (dataType == .dataTypes()$RAW())
            check(object@.data,  .requiredHMCols())
        else check(object@.data, .requiredDataCols())
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


#' @rdname dataType-methods
#' @aliases dataType,PerturbationData-method
#' @import data.table
setMethod(
    "dataType",
    signature = signature(obj="PerturbationData"),
    function(obj) obj@dataType)
