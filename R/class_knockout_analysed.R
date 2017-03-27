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


#' @title Data wrapper for analysed knockout data
#'
#' @name KnockoutAnalysis-class
#' @rdname knockout_analysis-class
#'
#' @description Abstract class \code{knockout.analysed} is a wrapper for a
#'   \code{data.table} object
#' containing the knockout data
#'
#' @slot .data the knockout data-set
#' @slot .inference the method for inferenced that has been used
setClass(
  "knockout.analysed",
  contains = "VIRTUAL",
  slots    = list(.data="data.table", .inference="character"),
  validity = function(object)
  {
    stopifnot(object@.inference %in% .inference.types())
  }
)

#' Data wrapper for analysed knockout data using an LMM
#'
#' @name LMMAnalysis-class
#' @rdname lmm_kockout_analysis-class
#'
#' @description Class \code{knockout.analysed.lmm} is a wrapper for a
#'   \code{data.table} object
#' containing the knockout data
#'
#' @slot .is.bootstrapped boolean whether bootstrapping has been done or not
setClass(
  "knockout.analysed.lmm",
  contains  = "knockout.analysed",
  slots     = list(.is.bootstrapped="logical"),
  prototype = prototype(.is.bootstrapped=FALSE,
                        .inference=.inference.types()$MIXED.MODEL)
)

#' Data wrapper for analysed knockout data using network diffusion
#'
#' @name DiffusionAnalysis-class
#' @rdname diffusion_kockout_analysis-class
#'
#' @description Class \code{knockout.analysed.diffusion} is a wrapper for a
#'   \code{data.table} object containing the knockout data
#'
#' @slot .is.bootstrapped boolean whether bootstrapping has been done or not
setClass(
  "knockout.analysed.diffusion",
  contains  = "knockout.analysed",
  slots     = list(.is.bootstrapped="logical"),
  prototype = prototype(.is.bootstrapped=FALSE)
)
