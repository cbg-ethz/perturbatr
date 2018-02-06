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


#' @title Smooth the results from an analysis using network diffusion
#'
#' @description TODO
#'
#' @export
#' @docType methods
#' @rdname diffuse-methods
#'
#' @import data.table
#'
#' @param obj  an analysed object
#' @param path   path to the network file (if \code{graph} is \code{NULL})
#' @param graph  an weighted adjacency matrix (if \code{path} is \code{NULL})
#' @param r  restart probability of the random if method \code{mrw} is selected
#' @param delete.nodes.on.degree  delete nodes from the graph with a degree of
#'  less or equal than \code{delete.nodes.on.degree}
#' @param do.bootstrap  run a diffusion on every bootstrap sample in case
#'  bootstrap samples are available
#' @param ...  additional parameters
#'
#' @return returns a \code{perturbation.diffusion.analysed} object
#'
#' @examples
#'  data(rnaiscreen)
#'  rnaiscreen.normalized <- preprocess(rnaiscreen, normalize="robust-z.score")
#'  res                   <- hm(rnaiscreen.normalized, effect.size=0.01)
#'  graph.file <- system.file("extdata", "graph_file.tsv", package = "perturbatr")
#'  diffu      <- diffuse(res, path=graph.file, r=0.1)
setGeneric(
  "diffuse",
  function(obj,
           path=NULL,
           graph=NULL,
           r=0.5,
           delete.nodes.on.degree=0,
           do.bootstrap=FALSE,
           ...)
  {
    standardGeneric("diffuse")
  },
  package="perturbation"
)

#' @rdname diffuse-methods
#' @aliases diffuse,perturbation.hm.analysed-method
#' @import data.table
#' @importFrom dplyr select filter
setMethod(
  "diffuse",
  signature=signature(obj="perturbation.hm.analysed"),
  function(obj,
           path=NULL,
           graph=NULL,
           r=0.5,
           delete.nodes.on.degree=0,
           do.bootstrap=FALSE,
           ...)
  {
    hits <- dplyr::select(obj@.gene.hits, GeneSymbol, Effect) %>%
      dplyr::mutate(Effect = abs(Effect))
    if (nrow(hits) == 0)
      stop("Your prior analysis did not yield hits for genes")

    bootstrap.hits <- NULL
    if (obj@.is.bootstrapped)
    {
      bootstrap.hits <- obj@.model.fit$ge.fdrs$ret %>%
        dplyr::filter(GeneSymbol %in% hits$GeneSymbol) %>%
        dplyr::select(-Mean, -Pval, -Qval, -Lower, -Upper)
    }

    ret   <- .diffuse(hits,
                      mod=obj,
                      bootstrap.hits=bootstrap.hits,
                      path=path,
                      graph=graph,
                      method=method,
                      r=r,
                      delete.nodes.on.degree=delete.nodes.on.degree,
                      do.bootstrap=do.bootstrap)
   ret
  }
)

#' @noRd
#' @import data.table igraph
.diffuse <- function(hits,
                     mod,
                     bootstrap.hits,
                     path,
                     graph,
                     method,
                     r,
                     delete.nodes.on.degree,
                     do.bootstrap)
{
  graph <- .get.graph(path, graph, delete.nodes.on.degree)
  adjm   <- igraph::get.adjacency(graph, attr="weight")

  l <- mrw(hits=hits,
           delete.nodes.on.degree=delete.nodes.on.degree,
           mod=mod,
           bootstrap.hits=bootstrap.hits,
           adjm=adjm,
           r=r,
           graph=graph,
           do.bootstrap=do.bootstrap)

  l
}

.get.graph <- function(path, graph, delete.nodes.on.degree)
{
  graph <- .read.graph(path=path, graph=graph)
  if (igraph::is.directed(graph))
    stop("Please provide an undirected graph")
  # get connected components
  comps <- igraph::components(graph)
  if (length(comps$csize) > 1)
    message("Only taking largest connected component to ensure ergodicity.")
  # get the genes that are not in the largest component
  non.max.comp.genes <- names(which(comps$membership != which.max(comps$csize)))
  # remove the genes that are not in the largest component
  # this is needed to ensure ergocity
  graph <- igraph::delete.vertices(graph, non.max.comp.genes)
  # delete vertexes with node degree less than ...
  graph <- igraph::delete.vertices(
    graph, igraph::V(graph)[igraph::degree(graph) <= delete.nodes.on.degree])

  graph
}
