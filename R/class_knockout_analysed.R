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


#' Make igraph recognizable by S4
setOldClass("igraph")


#' @title Data wrapper for analysed knockout data
#'
#' @name KnockoutAnalysis-class
#' @rdname knockout_analysis-class
#'
#' @description Abstract class \code{knockout.analysed} is a wrapper for a
#'   \code{data.table} object
#' containing the knockout data
#'
#' @slot .data  the knockout data-set
#' @slot .inference  the method for inferenced that has been used
#' @slot .is.bootstrapped  boolean whether bootstrap intervals have been
#' @slot .params  list of some used parameters
#'  created or not
setClass(
  "knockout.analysed",
  contains = "VIRTUAL",
  slots    = list(.data="data.table",
                  .params="list",
                  .inference="character",
                  .is.bootstrapped="logical"),
  validity = function(object)
  {
    stopifnot(object@.inference %in% .inference.types())
  }
)


#' Data wrapper for analysed knockout data using a standard hypothesis test
#'
#' @name TstatisticAnalysis-class
#' @rdname tstatisic_knockout_analysis-class
#'
#' @description Class \code{knockout.tstatistic.analysed} is a wrapper for a
#'   \code{data.table} object containing the knockout data
#'
#' @slot .gene.hits  prioritized genes
#'
setClass(
  "knockout.tstatistic.analysed",
  contains  = "knockout.analysed",
  slots     = list(.gene.hits="data.table"),
  prototype = prototype(.inference=.inference.types()$T.TEST,
                        .is.bootstrapped=FALSE)
)

#' Data wrapper for analysed knockout data using a standard hypothesis test
#'
#' @name HyperAnalysis-class
#' @rdname hyper_knockout_analysis-class
#'
#' @description Class \code{knockout.hyper.analysed} is a wrapper for a
#'   \code{data.table} object containing the knockout data
#'
#' @slot .gene.hits  prioritized genes
#'
setClass(
  "knockout.hyper.analysed",
  contains  = "knockout.analysed",
  slots     = list(.gene.hits="data.table"),
  prototype = prototype(.inference=.inference.types()$HYPERGEOMETRIC.TEST,
                        .is.bootstrapped=FALSE)
)


#' Data wrapper for analysed knockout data using an LMM
#'
#' @name LMMAnalysis-class
#' @rdname lmm_knockout_analysis-class
#'
#' @description Class \code{knockout.lmm.analysed} is a wrapper for a
#'   \code{data.table} object containing the knockout data
#'
#' @slot .gene.effects  the estimated effect sizes for genes
#' @slot .gene.pathogen.effects  the estimated effect sizes for genes on a
#'  viral level
#' @slot .infectivity.effects the  estimated effect sizes for different
#'  infectivity levels
#' @slot .gene.hits  prioritized genes
#' @slot .gene.pathogen.hits  prioritized genes on a viral level
#' @slot .model.fit  the fitted model with gene fdrs and gene-pathogen
#'  fdrs
setClass(
  "knockout.lmm.analysed",
  contains  = "knockout.analysed",
  slots     = list(.gene.effects          = "data.table",
                   .gene.pathogen.effects = "data.table",
                   .screen.type.effects   = "data.table",
                   .gene.hits             = "data.table",
                   .gene.pathogen.hits    = "data.table",
                   .model.fit             = "list"),
  prototype = prototype(.inference=.inference.types()$MIXED.MODEL)
)


#' @title Data wrapper for analysed knockout data using network diffusion
#'
#' @name DiffusionAnalysis-class
#' @rdname diffusion_knockout_analysis-class
#'
#' @import igraph
#'
#' @description Class \code{knockout.diffusion.analysed} is a wrapper for a
#'   \code{data.table} object containing the knockout data
#'
#' @slot .parameters  the parameters used for network analysis
#' @slot .intitial.model  the model that was provided for analysis
#' @slot .graph  an igraph object that served for the diffusion process
setClass(
  "knockout.diffusion.analysed",
  contains  = "knockout.analysed",
  slots     = list(.parameters    = "list",
                   .initial.model = "ANY",
                   .graph         = "igraph")
)
