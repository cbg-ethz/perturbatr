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
#' @slot geneEffects  the estimated effect sizes for genes
#' @slot inference  the method for inferenced that has been used
#' @slot isBootstrapped  boolean whether bootstrap intervals have been
#' @slot params  list of some used parameters
#'
setClass(
  "AbstractAnalysedPerturbationData",
  contains = "VIRTUAL",
  slots    = list(dataSet        = "tbl_df",
                  params         = "list",
                  geneEffects    = "tbl_df",
                  inference      = "character",
                  isBootstrapped = "logical"),
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
#'  using a hierarchical model. Class \code{HMAnalysedPerturbationData} exposes
#'  getters for its members of the same name, but no setters, because the data
#'  should be treated as
#'  constant once set. Objects of class \code{HMAnalysedPerturbationData}
#'  do not need to be
#'  constructed manually but are returned from calling \code{\link{hm}} (see
#'  the examples).
#'
#' @slot nestedGeneEffects  the estimated effect sizes for genes on a
#'  viral level
#' @slot modelFit  the fitted model
#'
#' @examples
#'  data(rnaiscreen)
#'  res <- hm(rnaiscreen, effect.size=0.01)
#'  class(res)
setClass(
    "HMAnalysedPerturbationData",
    contains  = c("AbstractAnalysedPerturbationData"),
    slots     = list(nestedGeneEffects = "tbl_df",
                     modelFit          = "list"),
    prototype = prototype(inference=inferenceTypes()$MIXED.MODEL)
)


#' @title Data wrapper for analysed perturbation data using network diffusion
#'
#' @name NetworkAnalysedPerturbationData-class
#' @rdname NetworkAnalysedPerturbationData-class
#' @exportClass NetworkAnalysedPerturbationData
#'
#' @description Class \code{NetworkAnalysedPerturbationData} is a wrapper for
#'  various objects of an analysis of a perturbation experiment done
#'  using network diffusion. Class \code{NetworkAnalysedPerturbationData}
#' exposes getters for its members of the same name,
#'  but no setters, because the data should be treated as
#'  constant once set. Objects of class \code{NetworkAnalysedPerturbationData}
#'  do not need to be
#'  constructed manually but are returned from calling \code{\link{diffuse}}
#'  (see the examples).
#'
#' @slot initialModel  the model that was provided for analysis
#' @slot graph  an igraph object that served for the diffusion process
#'
#' @examples
#'  data(rnaiscreen)
#'  hm.fit <- hm(rnaiscreen, effect.size=0.01)
#'  graph.file <- system.file("extdata", "graph_file.tsv",
#'                                        package = "perturbatr")
#'  res <- diffuse(hm.fit, path=graph.file, r=1)
setClass(
  "NetworkAnalysedPerturbationData",
  contains  = c("AbstractAnalysedPerturbationData"),
  slots     = list(initialModel = "ANY",
                   graph        = "igraph"),
  prototype = prototype(inference=inferenceTypes()$MRW.DIFFUSION,
                        isBootstrapped=FALSE)
)


#' @rdname dataSet-methods
#' @aliases dataSet,AbstractAnalysedPerturbationData-method
#' @import tibble
setMethod(
    "dataSet",
    signature = signature(obj="AbstractAnalysedPerturbationData"),
    function(obj) obj@dataSet)


#' @rdname params-methods
#' @aliases params,AbstractAnalysedPerturbationData-method
#' @import tibble
setMethod(
    "params",
    signature = signature(obj="AbstractAnalysedPerturbationData"),
    function(obj) obj@params)


#' @rdname inference-methods
#' @aliases inference,AbstractAnalysedPerturbationData-method
#' @import tibble
setMethod(
    "inference",
    signature = signature(obj="AbstractAnalysedPerturbationData"),
    function(obj) obj@inference)


#' @rdname isBootstrapped-methods
#' @aliases isBootstrapped,AbstractAnalysedPerturbationData-method
#' @import tibble
setMethod(
    "isBootstrapped",
    signature = signature(obj="AbstractAnalysedPerturbationData"),
    function(obj) obj@isBootstrapped)


#' @rdname geneEffects-methods
#' @aliases geneEffects,AbstractAnalysedPerturbationData-method
#' @import tibble
setMethod(
    "geneEffects",
    signature = signature(obj="AbstractAnalysedPerturbationData"),
    function(obj) obj@geneEffects)


#' @rdname nestedGeneEffects-methods
#' @aliases nestedGeneEffects,HMAnalysedPerturbationData-method
#' @import tibble
setMethod(
    "nestedGeneEffects",
    signature = signature(obj="HMAnalysedPerturbationData"),
    function(obj) obj@nestedGeneEffects)


#' @rdname modelFit-methods
#' @aliases modelFit,HMAnalysedPerturbationData-method
#' @import tibble
setMethod(
  "modelFit",
  signature = signature(obj="HMAnalysedPerturbationData"),
  function(obj) obj@modelFit)


#' @rdname graph-methods
#' @aliases graph,NetworkAnalysedPerturbationData-method
#' @import tibble
setMethod(
    "graph",
    signature = signature(obj="NetworkAnalysedPerturbationData"),
    function(obj) obj@graph)
