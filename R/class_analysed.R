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
#' @slot dataSet  the perturbation data-set
#' @slot inference  the method for inferenced that has been used
#' @slot isBootstrapped  boolean whether bootstrap intervals have been
#' @slot params  list of some used parameters
#'
setClass(
    "AbstractAnalysedPerturbationExperiment",
    contains = "VIRTUAL",
    slots    = list(dataSet="data.table",
                    params="list",
                    inference="character",
                    isBootstrapped="logical"),
    validity = function(object)
    {
        stopifnot(object@.inference %in% .inference.types())
    }
)


#' @rdname dataSet-methods
#' @aliases dataSet,AbstractAnalysedPerturbationExperiment-method
#' @import data.table
setMethod(
    "dataSet",
    signature = signature(obj="AbstractAnalysedPerturbationExperiment"),
    function(obj) obj@dataSet)


#' @rdname params-methods
#' @aliases params,AbstractAnalysedPerturbationExperiment-method
#' @import data.table
setMethod(
    "params",
    signature = signature(obj="AbstractAnalysedPerturbationExperiment"),
    function(obj) obj@params)


#' @rdname inference-methods
#' @aliases inference,AbstractAnalysedPerturbationExperiment-method
#' @import data.table
setMethod(
    "inference",
    signature = signature(obj="AbstractAnalysedPerturbationExperiment"),
    function(obj) obj@inference)


#' @rdname isBootstrapped-methods
#' @aliases isBootstrapped,AbstractAnalysedPerturbationExperiment-method
#' @import data.table
setMethod(
    "isBootstrapped",
    signature = signature(obj="AbstractAnalysedPerturbationExperiment"),
    function(obj) obj@isBootstrapped)


#' Data wrapper for analysed perturbation data using a standard test statistic
#'
#' @name AnalysedPerturbationExperiment-class
#' @rdname AnalysedPerturbationExperiment-class
#'
#' @description Class \code{AnalysedPerturbationExperiment} is a wrapper of an
#'  analysed data set using a test statistic such as a t-test or hypergeometric
#'  test.
#'
#' @slot geneHits  prioritized genes
#'
setClass(
    "AnalysedPerturbationExperiment",
    contains  = "AbstractAnalysedPerturbationExperiment",
    slots     = list(geneHits="data.table"),
    prototype = prototype(isBootstrapped=FALSE)
)


#' @rdname geneHits-methods
#' @aliases geneHits,AnalysedPerturbationExperiment-method
#' @import data.table
setMethod(
    "geneHits",
    signature = signature(obj="AnalysedPerturbationExperiment"),
    function(obj) obj@geneHits)


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
#' @slot geneEffects  the estimated effect sizes for genes
#' @slot nestedGeneEffects  the estimated effect sizes for genes on a
#'  viral level
#' @slot nestedGeneHits  nested prioritized genes
#' @slot modelFit  the fitted model
#'
setClass(
  "HMAnalysedPerturbationExperiment",
  contains  = "AnalysedPerturbationExperiment",
  slots     = list(geneEffects       = "data.table",
                   nestedGeneEffects = "data.table",
                   nestedGeneHits    = "data.table",
                   modelFit          = "list"),
  prototype = prototype(inference=.inference.types()$MIXED.MODEL)
)


#' @rdname geneEffects-methods
#' @aliases geneEffects,HMAnalysedPerturbationExperiment-method
#' @import data.table
setMethod(
    "geneEffects",
    signature = signature(obj="HMAnalysedPerturbationExperiment"),
    function(obj) obj@geneEffects)


#' @rdname nestedGeneEffects-methods
#' @aliases nestedGeneEffects,HMAnalysedPerturbationExperiment-method
#' @import data.table
setMethod(
    "nestedGeneEffects",
    signature = signature(obj="HMAnalysedPerturbationExperiment"),
    function(obj) obj@nestedGeneEffects)


#' @rdname nestedGeneHits-methods
#' @aliases nestedGeneHits,HMAnalysedPerturbationExperiment-method
#' @import data.table
setMethod(
    "nestedGeneHits",
    signature = signature(obj="HMAnalysedPerturbationExperiment"),
    function(obj) obj@nestedGeneHits)


#' @rdname modelFit-methods
#' @aliases modelFit,HMAnalysedPerturbationExperiment-method
#' @import data.table
setMethod(
    "modelFit",
    signature = signature(obj="HMAnalysedPerturbationExperiment"),
    function(obj) obj@modelFit)


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
    slots     = list(initialModel = "ANY",
                     graph        = "igraph"),
                     prototype = prototype(
                            inference=.inference.types()$MRW.DIFFUSION,
                            isBootstrapped=FALSE)
)


#' @rdname graph-methods
#' @aliases graph,NetworkAnalysedPerturbationExperiment-method
#' @import data.table
setMethod(
    "graph",
    signature = signature(obj="NetworkAnalysedPerturbationExperiment"),
    function(obj) obj@graph)
