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
#' @include class_data.R
#' @include class_analysed.R


#' @title Network diffusion
#'
#' @description Propagate the estimated gene effects from a previous analysis
#'  over a network using network diffusion. First the estimated effects are
#'  normalized and mapped to a given genetic network, for instance a PPI or
#'  co-expression network. Then the normalized effects are propagated across the
#'  edges of the network using a Markov random walk with restarts.
#'  By that the initial ranking of genes
#'  (as given by their absolute effect sizes) is re-evaluated and the genes are
#'  reordered. Thus network diffusion potentially reduced false negative hits.
#'
#' @export
#' @docType methods
#' @rdname diffuse-methods
#'
#' @import tibble
#'
#' @param obj  \code{HMAnalysedPerturbationData} object
#' @param graph a \code{data.frame} or \code{tibble} with three columns
#'  representing a symbolic edge list.
#'  The first two columns contain node ids. The third column has to be called
#'  \emph{weight} and is the \emph{non-negative} weight of the edge between
#'  the two nodes.
#' @param r  restart probability of the random walk
#' @param delete.nodes.on.degree  delete nodes from the graph with a degree of
#'  less or equal than \code{delete.nodes.on.degree}
#' @param do.bootstrap  run a diffusion on every bootstrap sample in case
#'  bootstrap samples are available
#' @param take.largest.component  if \code{true} takes only the largest
#'  connected component of the graph and discards all nodes that are not in
#'  the largest component. If \code{false} takes the compete graph.
#'
#' @return returns a \code{NetworkAnalysedPerturbationData} object
#'
#' @examples
#'  data(rnaiscreen)
#'  res        <- hm(rnaiscreen)
#'
#'  graph <- readRDS(system.file(
#'    "extdata", "graph_small.rds", package = "perturbatr"))
#'  res <- diffuse(hm.fit, graph=graph, r=1)
#'
setGeneric(
  "diffuse",
  function(obj,
           graph=NULL,
           r=0.5,
           delete.nodes.on.degree=0,
           do.bootstrap=FALSE,
           take.largest.component=TRUE)
  {
      standardGeneric("diffuse")
  }
)


#' @rdname diffuse-methods
#' @aliases diffuse,HMAnalysedPerturbationData-method
#' @import tibble
#' @importFrom dplyr select filter
#' @importFrom rlang .data
setMethod(
  "diffuse",
  signature=signature(obj="HMAnalysedPerturbationData"),
  function(obj,
           graph=NULL,
           r=0.5,
           delete.nodes.on.degree=0,
           do.bootstrap=FALSE,
           take.largest.component=TRUE)
  {
    hits <- dplyr::select(geneEffects(obj), .data$GeneSymbol, .data$Effect)
    hits <- dplyr::mutate(hits, "Effect" = abs(.data$Effect))
    if (nrow(hits) == 0)
        stop("Your prior analysis did not yield gene effect sizes")

    bootstrap.hits <- NULL
    if (isBootstrapped(obj))
    {
      bootstrap.hits <- modelFit(obj)$ge.fdrs$ret
      bootstrap.hits <- dplyr::filter(bootstrap.hits,
        .data$GeneSymbol %in% hits$GeneSymbol)
      bootstrap.hits <- dplyr::select(
        bootstrap.hits, -.data$Mean, -.data$Pval, -.data$Qval,
        -.data$Lower, -.data$Upper)
    }

    ret <- .diffuse(hits,
                    mod=obj,
                    bootstrap.hits=bootstrap.hits,
                    graph=graph,
                    r=r,
                    delete.nodes.on.degree=delete.nodes.on.degree,
                    do.bootstrap=do.bootstrap,
                    take.largest.component=take.largest.component)
    ret
  }
)


#' @noRd
#' @importFrom igraph get.adjacency
.diffuse <- function(hits,
                     mod,
                     bootstrap.hits,
                     graph,
                     r,
                     delete.nodes.on.degree,
                     do.bootstrap,
                     take.largest.component)
{
  graph <- .get.graph(graph, delete.nodes.on.degree, take.largest.component)
  adjm  <- igraph::get.adjacency(graph, attr="weight")

  l <- mrw(hits=hits,
           delete.nodes.on.degree=delete.nodes.on.degree,
           mod=mod,
           bootstrap.hits=bootstrap.hits,
           adjm=adjm,
           r=r,
           graph=graph,
           do.bootstrap=do.bootstrap,
           take.largest.component=take.largest.component)

  l
}


#' @noRd
#' @importFrom igraph is.directed components delete.vertices
#' @importFrom igraph V degree
.get.graph <- function(graph,delete.nodes.on.degree, take.largest.component)
{
  graph <- read.graph(graph)

  if (take.largest.component) {
    # get connected components
    comps <- igraph::components(graph)
    if (length(comps$csize) > 1) {
      message("Only taking largest connected component to ensure ergodicity.")
    }

    # get the genes that are not in the largest component
    non.max.comp.genes <- names(
      which(comps$membership != which.max(comps$csize)))
    # remove the genes that are not in the largest component
    # this is needed to ensure ergocity
    graph <- igraph::delete.vertices(graph, non.max.comp.genes)
  }

  # delete vertexes with node degree less than delete.nodes.on.degree
  graph <- igraph::delete.vertices(
      graph,
      igraph::V(graph)[igraph::degree(graph) <= delete.nodes.on.degree])

  graph
}
