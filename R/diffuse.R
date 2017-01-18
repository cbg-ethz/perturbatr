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


#' Extend the results from a <code>svd.prioritized</code> object by network diffusion.
#'
#' @export
#'
#' @import data.table
#'
#' @param obj  an \code{svd.prioritized} object
#' @param method  method that should be used for diffusion
#' \itemize{
#'   \item{neighbors }{ just looks at the neighbors :)}
#'   \item{mrw }{ do a Markov random walk}
#' }
#' @param path   path to the network file
#' @param r  restart probability of the random if method \code{mrw} is selected
#' @param node.start.count  number of nodes that are used to do the neighbors
#'  search if method \code{neighbors} is selected
#'  If the number of hits in \code{obj} exceeds \code{node.start.count},
#'  then the elements with the highest absolute effects are chosen. If this
#'  behaviour is not desired filter \code{obj} before
#' @param search.depth  how deep should the neighbor search go if method \code{neighbors} is selected
#' @param ...  additional parameters
diffuse <- function(obj, method=c("neighbors", "mrw"), path,
                    r=0.5, node.start.count=25, search.depth=5, ...)
{
  UseMethod("diffuse")
}

#' @noRd
#' @export
#' @import data.table igraph
diffuse.svd.prioritized.pmm <- function(obj, method=c("neighbors", "mrw"),
                                        path, r=0.5, node.start.count=25,
                                        search.depth=5, ...)
{
  if (!file.exists(path))
    stop(paste("Can't find: ", path, "!", sep=""))
  hits <- obj$gene.hits %>%
    dplyr::select(GeneSymbol, abs(Effect))
  res   <- .diffuse(hits, path,
                    match.arg(method),
                    r=0.5, node.start.count=25, search.depth=5)
  invisible(res)
}

#' @noRd
#' @import data.table
#' @importFrom igraph get.adjacency
.diffuse <- function(hits, path, method, r=0.5,
                     node.start.count=25, search.depth=5)
{
  graph <- .read.graph(path)
  adjm  <-  igraph::get.adjacency(graph, attr="weight")
  switch(method,
         "neighbors"= .knn(hits, adjm, node.start.count, search.depth, graph),
         "mrw"      = .mrw(hits, adjm, r, graph),
         stop("No suitable method found"))
}

#' @noRd
#' @import data.table
#' @importFrom diffusr random.walk
.mrw <- function(hits, adjm, r, graph)
{
  diffuse.data <- .init.starting.distribution(hits, adjm)
  mrw <- diffusr::random.walk(abs(diffuse.data$frame$Effect),
                              as.matrix(diffuse.data$adjm), r)
  diffuse.data$frame$DiffusionEffect <- mrw
  res <- diffuse.data$frame
  li <- list(diffusion=res, lmm.hits=hits, graph=graph)
  class(li) <- c("svd.diffused.mrw", "svd.diffused", class(li))
  return(li)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr left_join
.init.starting.distribution <- function(hits, adjm)
{
  res  <- dplyr::left_join(data.table::data.table(GeneSymbol=colnames(adjm)),
                           hits, by="GeneSymbol")
  data.table::setDT(res)[is.na(Effect), Effect := 0]
  return(list(frame=res, adjm=adjm))
}

#' @noRd
#' @import data.table igraph
#' @importFrom dplyr filter mutate
.knn <- function(hits, adjm, node.start.count, search.depth, graph)
{
  diffuse.data <- .init.starting.indexes(hits, adjm, node.start.count)
  neighs <- diffusr::nearest.neighbors(
    as.integer(dplyr::filter(diffuse.data$frame, Select==T)$Idx),
    as.matrix(diffuse.data$adjm), search.depth)
  res <- diffuse.data$frame
  names(neighs) <- filter(diffuse.data$frame,
                       Idx %in% as.integer(names(neighs)))$GeneSymbol
  li <- list(data=res, neighbors=knn, graph=graph)
  class(li) <- c("svd.diffused.knn", "svd.diffused", class(res))
  return(li)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr left_join mutate select
.init.starting.indexes <- function(hits, adjm, node.start.count)
{
  res  <- dplyr::left_join(data.table::data.table(GeneSymbol=colnames(adjm)),
                           hits)
  data.table::setDT(res)[is.na(Effect), Effect := 0]
  res <- res %>%
    .[order(-abs(Effect))] %>%
    dplyr::mutate(Idx = 1:.N) %>%
    dplyr::mutate(Select = (Idx <= node.start.count))
  data.table::setDT(res)[is.na(Select), Select := FALSE]
  return(list(frame=res, adjm=adjm))
}
