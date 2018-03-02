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
.required.data.cols  <- function()
{
    c("Condition", "Replicate", "GeneSymbol", "Perturbation", "Readout")
}


#' @noRd
.required.hm.cols   <- function()
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
#'
#' @slot data  the data set as a \code{data.table} object
#' @slot datatype  the type that describes your data best, for instance
#'  \code{raw} or \code{normalized}
#'
setClass(
    "PerturbationData",
    slots    = list(data="data.table", datatype="character")
    validity = function(object)
    {
        if (datatype == data.types$RAW()) check(object@.data, 
                                                .required.hm.cols)_)
        else check(object@.data, .required.data.cols())
        return(TRUE)
    }
)
