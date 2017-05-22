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
#' @include class_knockout_data.R
#' @include class_knockout_analysed.R

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
#' @param method  method that should be used for diffusion
#' \itemize{
#'   \item{knn }{ use a nearest neighbor approach}
#'   \item{mrw }{ do a Markov random walk with restarts}
#' }
#' @param path   path to the network file (if \code{graph} is \code{NULL})
#' @param graph  an weighted adjacency matrix (if \code{path} is \code{NULL})
#' @param r  restart probability of the random if method \code{mrw} is selected
#' @param node.start.count  number of nodes that are used to do the neighbors
#'  search if method \code{knn} is selected.
#'  If the number of hits in \code{obj} exceeds \code{node.start.count},
#'  then the elements with the highest absolute effects are chosen. If this
#'  behaviour is not desired filter \code{obj} before.
#' @param search.depth  how deep should the neighbor search go if method
#'  \code{nearest.neighbors} is selected
#' @param delete.nodes.on.degree  delete nodes from the graph with a degree of
#'  less or equal than \code{delete.nodes.on.degree}
#' @param ...  additional parameters
setGeneric(
  "diffuse",
  function(obj,
           method=c("knn", "mrw"),
           path=NULL,
           graph=NULL,
           r=0.5,
           node.start.count=25,
           search.depth=5,
           delete.nodes.on.degree=0,
           ...)
  {
    standardGeneric("diffuse")
  },
  package="knockout"
)

#' @rdname diffuse-methods
#' @aliases diffuse,knockout.lmm.analysed-method
#' @import data.table
#' @importFrom dplyr select filter
setMethod(
  "diffuse",
  signature=signature(obj="knockout.lmm.analysed"),
  function(obj,
           method=c("knn", "mrw"),
           path=NULL,
           graph=NULL,
           r=0.5,
           node.start.count=25,
           search.depth=5,
           delete.nodes.on.degree=0,
           ...)
  {
    method <- match.arg(method)
    hits   <- dplyr::select(obj@.gene.hits, GeneSymbol, abs(Effect))
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
                      node.start.count=node.start.count,
                      search.depth=search.depth,
                      delete.nodes.on.degree=delete.nodes.on.degree)

   ret
  }
)

#' @noRd
#' @import data.table igraph
.diffuse <- function(hits, mod, bootstrap.hits, path, graph, method, r,
                     node.start.count, search.depth, delete.nodes.on.degree)
{
  graph <- .read.graph(path=path, graph=graph)
  if (igraph::is.directed(graph))
    stop("Please provide an undirected graph")

  # get connected components
  comps        <- igraph::components(graph)
  # get the genes that are not in the largest component
  non.max.comp.genes <- names(which(comps$membership != which.max(comps$csize)))
  # remove the genes that are not in the largest component
  # this is needed to ensure ergocity
  graph <- igraph::delete.vertices(graph, non.max.comp.genes)

  graph <- igraph::delete.vertices(
    graph, igraph::V(graph)[
      igraph::degree(graph) <= delete.nodes.on.degree  ])
  adjm  <- igraph::get.adjacency(graph, attr="weight")

  l <- switch(method,
         "knn" = knn(hits, bootstrap.hits, adjm,
                     node.start.count, search.depth, graph),
         "mrw" = mrw(hits,
                     delete.nodes.on.degree=delete.nodes.on.degree,
                     mod=mod,
                     bootstrap.hits=bootstrap.hits,
                     adjm=adjm,
                     r=r,
                     graph=graph),
         stop("No suitable method found. Pick either from [knn/mrw]."))

  l
}
