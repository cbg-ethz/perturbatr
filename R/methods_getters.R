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


#' @title Getter for a data set
#'
#' @description Returns the data set that underlies an S4 wrapper class
#'  as \code{data.table}.
#'
#' @export
#' @docType methods
#' @rdname dataSet-methods
#'
#' @param obj  the object for which you want to extract the underlying data
#' @return  returns a \code{data.table}.
#' @examples
#'  data(rnaiscreen)
#'  dataSet(rnaiscreen)
setGeneric("dataSet", function(obj) standardGeneric("dataSet"))


#' @title Getter for parameters used for analysis of perturbation data
#'
#' @description Returns the parameters used in the analysis of a perturbation
#'  analysis. Parameters are for instance the significance level, the effect
#'  size or restart probability of a random walk.
#'
#' @export
#' @docType methods
#' @rdname params-methods
#'
#' @param obj  the object for which you want to extract the underlying params
#' @return  returns a \code{list}.
#' @examples
#'  data(rnaiscreen)
#'  ft <- hm(rnaiscreen)
#'  params(ft)
setGeneric("params", function(obj) standardGeneric("params"))


#' @title Getter for inference used for analysis of perturbation data
#'
#' @description Returns the inference used in the analysis of a perturbation
#'  analysis. This can for instance be a standard t-test or a hierarchical
#'  model.
#'
#' @export
#' @docType methods
#' @rdname inference-methods
#'
#' @param obj  the object for which you want to extract the underlying inference
#' @return  returns a \code{character}.
#' @examples
#'  data(rnaiscreen)
#'  ft <- hm(rnaiscreen)
#'  inference(ft)
setGeneric("inference", function(obj) standardGeneric("inference"))


#' @title Getter for boolean if boostrapping was used
#'
#' @description Returns a boolean if for the analysed object bootstrapping was
#'  used to create confidence intervals.
#'
#' @export
#' @docType methods
#' @rdname isBootstrapped-methods
#'
#' @param obj  the object for which you want to extract the boolean
#' @return  returns a \code{boolean}.
#' @examples
#'  data(rnaiscreen)
#'  ft <- hm(rnaiscreen)
#'  isBootstrapped(ft)
setGeneric("isBootstrapped", function(obj) standardGeneric("isBootstrapped"))


#' @title Getter for identified essential genes
#'
#' @description Returns a \code{data.table} containing the genes that have been
#'  identified as essential genes in a perturbation screen.
#'
#' @export
#' @docType methods
#' @rdname geneHits-methods
#'
#' @param obj  the object for which you want to extract the underlying gene hits
#' @return  returns a \code{data.table}.
#' @examples
#'  data(rnaiscreen)
#'  ft <- hm(rnaiscreen)
#'  geneHits(ft)
setGeneric("geneHits", function(obj) standardGeneric("geneHits"))


#' @title Getter for the complete list of genes and their effect sizes
#'
#' @description Returns a \code{data.table} containing the list of genes that
#' have been used in the study along with their estimated effect sizes.
#'
#' @export
#' @docType methods
#' @rdname geneEffects-methods
#'
#' @param obj  the object for which you want to extract the underlying gene
#' @return  returns a \code{list}.
#' @examples
#'  data(rnaiscreen)
#'  ft <- hm(rnaiscreen)
#'  geneEffects(ft)
setGeneric("geneEffects", function(obj) standardGeneric("geneEffects"))


#' @title Getter for the completed list of nested gene effects
#'
#' @description Returns a \code{data.table} containing the list of nested genes
#'  effects, i.e. a table of nested cluster effects.
#'
#' @export
#' @docType methods
#' @rdname nestedGeneEffects-methods
#'
#' @param obj  the object for which you want to extract the underlying effects
#' @return  returns a \code{data.table}.
#' @examples
#'  data(rnaiscreen)
#'  ft <- hm(rnaiscreen)
#'  nestedGeneEffects(ft)
setGeneric(
    "nestedGeneEffects", function(obj) standardGeneric("nestedGeneEffects"))


#' @title Getter for identified nested essential genes
#'
#' @description Returns a \code{data.table} containing the genes that have been
#'  identified as essential genes in a perturbation screen.
#'
#' @export
#' @docType methods
#' @rdname nestedGeneHits-methods
#'
#' @param obj  the object for which you want to extract the underlying gene hits
#' @return  returns a \code{data.table}.
#' @examples
#'  data(rnaiscreen)
#'  ft <- hm(rnaiscreen)
#'  nestedGeneHits(ft)
setGeneric("nestedGeneHits", function(obj) standardGeneric("nestedGeneHits"))


#' @title Getter for model fit
#'
#' @description Returns the mixed effects model fit from the analysis using a
#'  hierarchical model.
#'
#' @export
#' @docType methods
#' @rdname modelFit-methods
#'
#' @param obj  the object for which you want to extract the underlying fit
#' @return  returns a \code{list}
#' @examples
#'  data(rnaiscreen)
#'  ft <- hm(rnaiscreen)
#'  modelFit(ft)
setGeneric("modelFit", function(obj) standardGeneric("modelFit"))


#' @title Getter for graph used for network diffusion
#'
#' @description Returns the graph that has been used for analysis using network
#'  diffusion.
#'
#' @export
#' @docType methods
#' @rdname graph-methods
#'
#' @param obj  the object for which you want to extract the underlying graph
#' @return  returns a \code{igraph} object
#' @examples
#'  \dontrun{
#'   data(rnaiscreen)
#'   graph.file <- system.file("extdata", "graph_file.tsv",
#'                             package = "perturbatr")
#'   ft <- hm(rnaiscreen, effect.size=0.01)
#'   diffu <- diffuse(ft, path=graph.file, r=0.1)
#'   graph(diffu)
#'  }
setGeneric("graph", function(obj) standardGeneric("graph"))
