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


#' Make igraph recognizable by S4
setOldClass("igraph")


#' @title Data wrapper for analysed perturbation data
#'
#' @name PerturbationAnalysis-class
#' @rdname perturbation_analysis-class
#'
#' @description Abstract class \code{perturbation.analysed} is a wrapper for a
#'   \code{data.table} object
#' containing the perturbation data
#'
#' @slot .data  the perturbation data-set
#' @slot .inference  the method for inferenced that has been used
#' @slot .is.bootstrapped  boolean whether bootstrap intervals have been
#' @slot .params  list of some used parameters
#'  created or not
setClass(
  "perturbation.analysed",
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


#' Data wrapper for analysed perturbation data using a hierarchical model
#'
#' @name HMAnalysis-class
#' @rdname hm_perturbation_analysis-class
#'
#' @description Class \code{perturbation.hm.analysed} is a wrapper for a
#'   \code{data.table} object containing the perturbation data
#'
#' @slot .gene.effects  the estimated effect sizes for genes
#' @slot .nested.gene.effects  the estimated effect sizes for genes on a
#'  viral level
#' @slot .gene.hits  prioritized genes
#' @slot .neted.gene.hits  nested prioritized genes
#' @slot .model.fit  the fitted model
setClass(
  "perturbation.hm.analysed",
  contains  = "perturbation.analysed",
  slots     = list(.gene.effects          = "data.table",
                   .nested.gene.effects   = "data.table",
                   .gene.hits             = "data.table",
                   .nested.gene.hits      = "data.table",
                   .model.fit             = "list"),
  prototype = prototype(.inference=.inference.types()$MIXED.MODEL)
)


#' @title Data wrapper for analysed perturbation data using network diffusion
#'
#' @name DiffusionAnalysis-class
#' @rdname diffusion_perturbation_analysis-class
#'
#' @import igraph
#'
#' @description Class \code{perturbation.diffusion.analysed} is a wrapper for a
#'   \code{data.table} object containing the perturbation data
#'
#' @slot .intitial.model  the model that was provided for analysis
#' @slot .graph  an igraph object that served for the diffusion process
setClass(
  "perturbation.diffusion.analysed",
  contains  = c("perturbation.analysed", "VIRTUAL"),
  slots     = list(.initial.model = "ANY",
                   .graph         = "igraph")
)



#' @title Data wrapper for analysed perturbation data using Markov random walks
#'  with restarts
#'
#' @name MRW-DiffusionAnalysis-class
#' @rdname mrw_diffusion_perturbation_analysis-class
#'
#' @import igraph
#'
#' @description Class \code{perturbation.mrw.diffusion.analysed} is a wrapper
#' for a \code{data.table} object containing the perturbation data
#'
setClass(
  "perturbation.mrw.diffusion.analysed",
  contains  = "perturbation.diffusion.analysed",
  prototype = prototype(.inference=.inference.types()$MRW.DIFFUSION,
                        .is.bootstrapped=FALSE)
)


#' Data wrapper for analysed perturbation data using a standard hypothesis test
#'
#' @name TstatisticAnalysis-class
#' @rdname tstatisic_perturbation_analysis-class
#'
#' @description Class \code{perturbation.tstatistic.analysed} is a wrapper for a
#'   \code{data.table} object containing the perturbation data
#'
#' @slot .gene.hits  prioritized genes
#'
setClass(
	"perturbation.tstatistic.analysed",
	contains  = "perturbation.analysed",
	slots     = list(.gene.hits="data.table"),
	prototype = prototype(.inference=.inference.types()$T.TEST,
												.is.bootstrapped=FALSE)
)


#' Data wrapper for analysed perturbation data using a standard hypothesis test
#'
#' @name HyperAnalysis-class
#' @rdname hyper_perturbation_analysis-class
#'
#' @description Class \code{perturbation.hyper.analysed} is a wrapper for a
#'   \code{data.table} object containing the perturbation data
#'
#' @slot .gene.hits  prioritized genes
#'
setClass(
	"perturbation.hyper.analysed",
	contains  = "perturbation.analysed",
	slots     = list(.gene.hits="data.table"),
	prototype = prototype(.inference=.inference.types()$HYPERGEOMETRIC.TEST,
												.is.bootstrapped=FALSE)
)
