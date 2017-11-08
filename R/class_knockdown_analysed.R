# knockdown: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2015 - 2016 Simon Dirmeier
#
# This file is part of knockdown
#
# knockdown is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# knockdown is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with knockdown. If not, see <http://www.gnu.org/licenses/>.


#' @include util_enums.R


#' Make igraph recognizable by S4
setOldClass("igraph")


#' @title Data wrapper for analysed knockdown data
#'
#' @name KnockdownAnalysis-class
#' @rdname knockdown_analysis-class
#'
#' @description Abstract class \code{knockdown.analysed} is a wrapper for a
#'   \code{data.table} object
#' containing the knockdown data
#'
#' @slot .data  the knockdown data-set
#' @slot .inference  the method for inferenced that has been used
#' @slot .is.bootstrapped  boolean whether bootstrap intervals have been
#' @slot .params  list of some used parameters
#'  created or not
setClass(
  "knockdown.analysed",
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


#' Data wrapper for analysed knockdown data using a standard hypothesis test
#'
#' @name TstatisticAnalysis-class
#' @rdname tstatisic_knockdown_analysis-class
#'
#' @description Class \code{knockdown.tstatistic.analysed} is a wrapper for a
#'   \code{data.table} object containing the knockdown data
#'
#' @slot .gene.hits  prioritized genes
#'
setClass(
  "knockdown.tstatistic.analysed",
  contains  = "knockdown.analysed",
  slots     = list(.gene.hits="data.table"),
  prototype = prototype(.inference=.inference.types()$T.TEST,
                        .is.bootstrapped=FALSE)
)

#' Data wrapper for analysed knockdown data using a standard hypothesis test
#'
#' @name ChisqStatisticAnalysis-class
#' @rdname chisq_statisic_knockdown_analysis-class
#'
#' @description Class \code{knockdown.chisqstatistic.analysed} is a wrapper for a
#'   \code{data.table} object containing the knockdown data
#'
#' @slot .gene.hits  prioritized genes
#'
setClass(
  "knockdown.chisqstatistic.analysed",
  contains  = "knockdown.analysed",
  slots     = list(.gene.hits="data.table"),
  prototype = prototype(.inference=.inference.types()$CHISQ.TEST,
                        .is.bootstrapped=FALSE)
)


#' Data wrapper for analysed knockdown data using a standard hypothesis test
#'
#' @name HyperAnalysis-class
#' @rdname hyper_knockdown_analysis-class
#'
#' @description Class \code{knockdown.hyper.analysed} is a wrapper for a
#'   \code{data.table} object containing the knockdown data
#'
#' @slot .gene.hits  prioritized genes
#'
setClass(
  "knockdown.hyper.analysed",
  contains  = "knockdown.analysed",
  slots     = list(.gene.hits="data.table"),
  prototype = prototype(.inference=.inference.types()$HYPERGEOMETRIC.TEST,
                        .is.bootstrapped=FALSE)
)


#' Data wrapper for analysed knockdown data using an LMM
#'
#' @name LMMAnalysis-class
#' @rdname lmm_knockdown_analysis-class
#'
#' @description Class \code{knockdown.lmm.analysed} is a wrapper for a
#'   \code{data.table} object containing the knockdown data
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
  "knockdown.lmm.analysed",
  contains  = "knockdown.analysed",
  slots     = list(.gene.effects          = "data.table",
                   .gene.pathogen.effects = "data.table",
                   .screen.type.effects   = "data.table",
                   .gene.hits             = "data.table",
                   .gene.pathogen.hits    = "data.table",
                   .model.fit             = "list"),
  prototype = prototype(.inference=.inference.types()$MIXED.MODEL)
)


#' @title Data wrapper for analysed knockdown data using network diffusion
#'
#' @name DiffusionAnalysis-class
#' @rdname diffusion_knockdown_analysis-class
#'
#' @import igraph
#'
#' @description Class \code{knockdown.diffusion.analysed} is a wrapper for a
#'   \code{data.table} object containing the knockdown data
#'
#' @slot .intitial.model  the model that was provided for analysis
#' @slot .graph  an igraph object that served for the diffusion process
setClass(
  "knockdown.diffusion.analysed",
  contains  = c("knockdown.analysed", "VIRTUAL"),
  slots     = list(.initial.model = "ANY",
                   .graph         = "igraph")
)



#' @title Data wrapper for analysed knockdown data using Markov random walks
#'  with restarts
#'
#' @name MRW-DiffusionAnalysis-class
#' @rdname mrw_diffusion_knockdown_analysis-class
#'
#' @import igraph
#'
#' @description Class \code{knockdown.mrw.diffusion.analysed} is a wrapper for a
#'   \code{data.table} object containing the knockdown data
#'
setClass(
  "knockdown.mrw.diffusion.analysed",
  contains  = "knockdown.diffusion.analysed",
  prototype = prototype(.inference.types()$MRW.DIFFUSION)
)


#' @title Data wrapper for analysed knockdown data using nearest neighbors
#'
#' @name kNN-DiffusionAnalysis-class
#' @rdname kNN_diffusion_knockdown_analysis-class
#'
#' @import igraph
#'
#' @description Class \code{knockdown.mrw.diffusion.analysed} is a wrapper for a
#'   \code{data.table} object containing the knockdown data
#'
#' @slot .neighbors  a \code{data.table} with found neighbors for every start
#'  gene
setClass(
  "knockdown.knn.diffusion.analysed",
  contains  = "knockdown.diffusion.analysed",
  slots     = list(.neighbors   = "data.table"),
  prototype = prototype(.inference.types()$NN.DIFFUSION,
                        .is.bootstrapped=FALSE)
)
