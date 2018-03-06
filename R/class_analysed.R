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
#' @include methods_getters.R


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
    "AbstractAnalysedPerturbationData",
    contains = "VIRTUAL",
    slots    = list(dataSet="data.table",
                    params="list",
                    inference="character",
                    isBootstrapped="logical"),
    validity = function(object)
    {
        stopifnot(object@inference %in% inferenceTypes())
    }
)


#' Data wrapper for analysed perturbation data using a hierarchical model
#'
#' @name HMAnalysedPerturbationData-class
#' @rdname HMAnalysedPerturbationData-class
#' @exportClass HMAnalysedPerturbationData
#'
#' @description Class \code{HMAnalysedPerturbationData} is a wrapper for
#'  various objects of an analysis of a perturbation experiment done
#'  using a hierarchical model.
#'
#' @slot geneEffects  the estimated effect sizes for genes
#' @slot geneHits  prioritized genes
#' @slot nestedGeneEffects  the estimated effect sizes for genes on a
#'  viral level
#' @slot nestedGeneHits  nested prioritized genes
#' @slot modelFit  the fitted model
#'
setClass(
  "HMAnalysedPerturbationData",
  contains  = c("AbstractAnalysedPerturbationData"),
  slots     = list(geneEffects       = "data.table",
                   geneHits          = "data.table",
                   nestedGeneEffects = "data.table",
                   nestedGeneHits    = "data.table",
                   modelFit          = "list"),
  prototype = prototype(inference=inferenceTypes()$MIXED.MODEL)
)


#' @title Data wrapper for analysed perturbation data using network diffusion
#'
#' @name NetworkAnalysedPerturbationData-class
#' @rdname NetworkAnalysedPerturbationData-class
#' @exportClass NetworkAnalysedPerturbationData
#'
#' @import igraph
#'
#' @description Class \code{NetworkAnalysedPerturbationData} is a wrapper for
#'  various objects of an analysis of a perturbation experiment done
#'  using network diffusion. See also \code{\link{diffuse}}.
#'
#' @slot intitialModel  the model that was provided for analysis
#' @slot graph  an igraph object that served for the diffusion process
setClass(
    "NetworkAnalysedPerturbationData",
    contains  = c("AbstractAnalysedPerturbationData", "VIRTUAL"),
    slots     = list(initialModel = "ANY",
                     graph        = "igraph"),
                     prototype = prototype(
                            inference=inferenceTypes()$MRW.DIFFUSION,
                            isBootstrapped=FALSE)
)


#' @rdname dataSet-methods
#' @aliases dataSet,AbstractAnalysedPerturbationData-method
#' @import data.table
setMethod(
    "dataSet",
    signature = signature(obj="AbstractAnalysedPerturbationData"),
    function(obj) obj@dataSet)


#' @rdname params-methods
#' @aliases params,AbstractAnalysedPerturbationData-method
#' @import data.table
setMethod(
    "params",
    signature = signature(obj="AbstractAnalysedPerturbationData"),
    function(obj) obj@params)


#' @rdname inference-methods
#' @aliases inference,AbstractAnalysedPerturbationData-method
#' @import data.table
setMethod(
    "inference",
    signature = signature(obj="AbstractAnalysedPerturbationData"),
    function(obj) obj@inference)


#' @rdname isBootstrapped-methods
#' @aliases isBootstrapped,AbstractAnalysedPerturbationData-method
#' @import data.table
setMethod(
    "isBootstrapped",
    signature = signature(obj="AbstractAnalysedPerturbationData"),
    function(obj) obj@isBootstrapped)


#' @rdname geneHits-methods
#' @aliases geneHits,HMAnalysedPerturbationData-method
#' @import data.table
setMethod(
    "geneHits",
    signature = signature(obj="HMAnalysedPerturbationData"),
    function(obj) obj@geneHits)


#' @rdname geneEffects-methods
#' @aliases geneEffects,HMAnalysedPerturbationData-method
#' @import data.table
setMethod(
    "geneEffects",
    signature = signature(obj="HMAnalysedPerturbationData"),
    function(obj) obj@geneEffects)


#' @rdname nestedGeneEffects-methods
#' @aliases nestedGeneEffects,HMAnalysedPerturbationData-method
#' @import data.table
setMethod(
    "nestedGeneEffects",
    signature = signature(obj="HMAnalysedPerturbationData"),
    function(obj) obj@nestedGeneEffects)


#' @rdname nestedGeneHits-methods
#' @aliases nestedGeneHits,HMAnalysedPerturbationData-method
#' @import data.table
setMethod(
    "nestedGeneHits",
    signature = signature(obj="HMAnalysedPerturbationData"),
    function(obj) obj@nestedGeneHits)


#' @rdname modelFit-methods
#' @aliases modelFit,HMAnalysedPerturbationData-method
#' @import data.table
setMethod(
  "modelFit",
  signature = signature(obj="HMAnalysedPerturbationData"),
  function(obj) obj@modelFit)


#' @rdname graph-methods
#' @aliases graph,NetworkAnalysedPerturbationData-method
#' @import data.table
setMethod(
    "graph",
    signature = signature(obj="NetworkAnalysedPerturbationData"),
    function(obj) obj@graph)
