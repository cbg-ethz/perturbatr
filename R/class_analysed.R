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
#' @noRd
#'
#' @slot data  the perturbation data-set
#' @slot inference  the method for inferenced that has been used
#' @slot is.bootstrapped  boolean whether bootstrap intervals have been
#' @slot params  list of some used parameters
#'
setClass(
    "AbstractAnalysedPerturbationExperiment",
    contains = "VIRTUAL",
    slots    = list(data="data.table",
                    params="list",
                    inference="character",
                    is.bootstrapped="logical"),
    validity = function(object)
    {
        stopifnot(object@.inference %in% .inference.types())
    }
)


#' Data wrapper for analysed perturbation data using a standard test statistic
#'
#' @name AnalysedPerturbationExperiment-class
#' @rdname AnalysedPerturbationExperiment-class
#'
#' @description Class \code{AnalysedPerturbationExperiment} is a wrapper of an
#'  analysed data set using a test statistic such as a t-test or hypergeometric
#'  test.
#'
#' @slot gene.hits  prioritized genes
#'
setClass(
    "AnalysedPerturbationExperiment",
    contains  = "AbstractAnalysedPerturbationExperiment",
    slots     = list(gene.hits="data.table"),
    prototype = prototype(is.bootstrapped=FALSE)
)


#' Data wrapper for analysed perturbation data using a hierarchical model
#'
#' @name HMAnalysedPerturbationExperiment-class
#' @rdname HMAnalysedPerturbationExperiment-class
#' @exportClass HMAnalysedPerturbationExperiment
#'
#' @description Class \code{pHMAnalysedPerturbationExperiment} is a wrapper for
#'  various objects of an analysis of a perturbation experiment done
#'  using a hierarchical model.
#'
#' @slot gene.effects  the estimated effect sizes for genes
#' @slot nested.gene.effects  the estimated effect sizes for genes on a
#'  viral level
#' @slot nested.gene.hits  nested prioritized genes
#' @slot model.fit  the fitted model
#'
setClass(
  "HMAnalysedPerturbationExperiment",
  contains  = "AnalysedPerturbationExperiment",
  slots     = list(gene.effects          = "data.table",
                   nested.gene.effects   = "data.table",
                   nested.gene.hits      = "data.table",
                   model.fit             = "list"),
  prototype = prototype(.inference=.inference.types()$MIXED.MODEL)
)


#' @title Data wrapper for analysed perturbation data using network diffusion in
#'
#' @name NetworkAnalysedPerturbationExperiment-class
#' @rdname NetworkAnalysedPerturbationExperiment-class
#'
#' @import igraph
#'
#' @description Class \code{DiffusionAnalysedPerturbationExperiment} is a
#'  various objects of an analysis of a perturbation experiment done
#'  using network diffusion. See also \code{\link{diffusion}}.
#'
#' @slot intitial.model  the model that was provided for analysis
#' @slot graph  an igraph object that served for the diffusion process
setClass(
  "NetworkAnalysedPerturbationExperiment",
  contains  = c("AbstractAnalysedPerturbationExperiment", "VIRTUAL"),
  slots     = list(initial.model = "ANY",
                   graph         = "igraph",
                   inference=.inference.types()$MRW.DIFFUSION,
                   is.bootstrapped=FALSE)
)
