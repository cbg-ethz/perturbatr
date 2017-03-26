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


#' Data wrapper for analysed knockout data
#'
#' @name AnalysedKnockoutData-class
#'
#' @description Abstract class \code{knockout.analysed} is a wrapper for a
#'   \code{data.table} object
#' containing the knockout data
#'
#' @slot .data the knockout data-set
#' @slot .inference the method for inferenced that has been used
setClass(
  "AnalysedKnockoutData",
  contains="VIRTUAL",
  slots     = list(.data="data.table", .inference="character"),
  validity=function(object)
  {
    stopifnot(object@.inference %in% .inference.types())
  }
)

#' Data wrapper for analysed knockout data using an LMM
#'
#' @name LMMAnalysedKnockoutData-class
#'
#' @description Class \code{knockout.analysed.lmm} is a wrapper for a
#'   \code{data.table} object
#' containing the knockout data
#'
#' @slot .is.bootstrapped boolean whether bootstrapping has been done or not
setClass(
  "LMMAnalysedKnockoutData",
  contains  = "AnalysedKnockoutData",
  slots     = list(.is.bootstrapped="logical"),
  prototype = prototype(.is.bootstrapped=FALSE,
                        .inference=.inference.types()$MIXED.MODEL)
)

#' Data wrapper for analysed knockout data using network diffusion
#'
#' @name DiffusionAnalysedKnockoutData-class
#'
#' @description Class \code{knockout.analysed.diffusion} is a wrapper for a
#'   \code{data.table} object
#' containing the knockout data
#'
#' @slot .is.bootstrapped boolean whether bootstrapping has been done or not
setClass(
  "DiffusionAnalysedKnockoutData",
  contains  = "AnalysedKnockoutData",
  slots     = list(.is.bootstrapped="logical"),
  prototype = prototype(.is.bootstrapped=FALSE)
)
