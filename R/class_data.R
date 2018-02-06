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
  c("Condition", "Replicate", "GeneSymbol", "Perturbation","Readout")
}

#' @noRd
.required.hm.cols   <- function()
{
  c("Condition", "GeneSymbol", "Weight", "Readout")
}

#' @noRd
#' @slot .data the perturbation data-set
setClass(
  "perturbation.data",
  contains = "VIRTUAL",
  slots    = list(.data="data.table")
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
  	.check(object@.data, .required.data.cols())
  	return(TRUE)
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
  	.check(object@.data, .required.data.cols())
    return(TRUE)
  }
)


#' @title Data wrapper for hierarchical model perturbation data
#'
#' @description Class \code{perturbation.hm.data} is a wrapper for a
#' \code{data.table} object containing the normalized data-set
#'
#' @name HM-data
#' @rdname perturbation_hm_data-class
#'
setClass(
  "perturbation.hm.data",
  contains  = "perturbation.data",
  validity  = function(object)
  {
  	.check(object@.data, .required.hm.cols())
    return(TRUE)
  }
)
