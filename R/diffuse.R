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
#' @param delete.nodes.on.degree  delete nodes from the graph with a degree of less that \code{delete.nodes.on.degree}
#' @param ...  additional parameters
diffuse <- function(obj, method=c("neighbors", "mrw"), path,
                    r=0.5, node.start.count=25, search.depth=5,
                    delete.nodes.on.degree, ...)
{
  UseMethod("diffuse")
}

#' @noRd
#' @export
#' @import data.table igraph
#' @importFrom dplyr filter select
#' @method diffuse svd.prioritized.pmm
diffuse.svd.prioritized.pmm <- function(obj, method=c("neighbors", "mrw"),
                                        path, r=0.5, node.start.count=25,
                                        search.depth=5, delete.nodes.on.degree=1,
                                        ...)
{
  if (!file.exists(path))
    stop(paste("Can't find: ", path, "!", sep=""))
  hits <- obj$gene.hits %>%
    dplyr::select(GeneSymbol, abs(Effect))

  # TODO: optimize
  loocv.hits <- obj$fit$fit$gene.fdrs %>%
    dplyr::filter(GeneSymbol %in% hits$GeneSymbol) %>%
    dplyr::select(-MeanLOOCVEffect, -Pval, -FDR)

  res   <- .diffuse(hits, loocv.hits,
                    path, match.arg(method),
                    r=r, node.start.count=node.start.count, search.depth=search.depth,
                    delete.nodes.on.degree=delete.nodes.on.degree)
  invisible(res)
}

#' @export
#' @import data.table igraph
#' @importFrom dplyr filter select rename
#' @method diffuse data.table
diffuse.data.table <- function(obj, method=c("neighbors", "mrw"),
                               path, r=0.5, node.start.count=25,
                               search.depth=5, delete.nodes.on.degree=1,
                               ...)
{
  if (!file.exists(path))
    stop(paste("Can't find: ", path, "!", sep=""))
  if ("Readout" %in% colnames(obj))
  {
    hits <- obj %>%
      dplyr::select(GeneSymbol, abs(Readout)) %>%
      dplyr::rename(Effect=Readout)
  }
  else if ("MeanEffect" %in% colnames(obj))
  {
    hits <- obj %>%
      dplyr::select(GeneSymbol, abs(MeanEffect)) %>%
      dplyr::rename(Effect=MeanEffect)
  }
  else if ("Effect" %in% colnames(obj))
  {
    hits <- obj %>% dplyr::select(GeneSymbol, abs(Effect))
  }
  else
  {
    stop("Neither 'Readout' nor 'Effect' found in colnames")
  }
  res   <- .diffuse(hits=hits,
                    loocv.hits=NULL,
                    path=path,
                    method=match.arg(method),
                    r=r,
                    node.start.count=node.start.count,
                    search.depth=search.depth,
                    delete.nodes.on.degree=delete.nodes.on.degree)
  invisible(res)
}

#' @noRd
#' @import data.table igraph
.diffuse <- function(hits, loocv.hits, path, method, r,
                     node.start.count, search.depth, delete.nodes.on.degree)
{
  graph <- .read.graph(path)
  graph <- igraph::delete.vertices(
    graph, igraph::V(graph)[
      igraph::degree(graph) <= delete.nodes.on.degree  ])
  adjm  <- igraph::get.adjacency(graph, attr="weight")
  switch(method,
         "neighbors"= .knn(hits, loocv.hits, adjm, node.start.count, search.depth, graph),
         "mrw"      = .mrw(hits, loocv.hits, adjm, r, graph),
         stop("No suitable method found"))
}

#' @noRd
#' @import data.table
#' @importFrom diffusr random.walk
.mrw <- function(hits, loocv.hits, adjm, r, graph)
{
  diffuse.data <- .do.mrw(hits, adjm, r)
  if (!is.null(loocv.hits)) {
    pvals <- .significance.mrw(loocv.hits, adjm, r)
    res   <- dplyr::left_join(diffuse.data$frame, pvals, by="GeneSymbol")
  }
  else {
    res  <- diffuse.data$frame
  }
  li  <- list(diffusion=dplyr::select(res, GeneSymbol, Effect, DiffusionEffect),
              diffusion.model=res,
              lmm.hits=hits,
              graph=graph)
  class(li) <- c("svd.diffused.mrw", "svd.diffused")
  return(li)
}

.do.mrw <- function(hits, adjm, r)
{
  diffuse.data <- .init.starting.distribution(hits, adjm)
  mrw <- diffusr::random.walk(abs(diffuse.data$frame$Effect),
                              as.matrix(diffuse.data$adjm), r)
  diffuse.data$frame$DiffusionEffect <- mrw
  diffuse.data
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
#' @importFrom tidyr gather
.significance.mrw <- function(loocv.hits, adjm, r)
{
  loocv.g <- tidyr::gather(loocv.hits, LOOCV, Effect, -GeneSymbol) %>%
    as.data.table
  li <- list()
  # TODO parallalize
  for (lo in unique(loocv.g$LOOCV))
  {
    hits <- dplyr::filter(loocv.g, LOOCV==lo) %>% dplyr::select(-LOOCV)
    dd.lo <- .do.mrw(hits, adjm, r)
    dd.lo$frame$loocv <- lo
    li[[lo]] <- dd.lo$frame
  }
  # TODO: how to do hypothesis test here?
  # -> calculate means -> calculate alphas -> calculate variance -> calculate alpha_0
  flat.dat <- do.call("rbind", lapply(li, function(e) e)) %>% dplyr::select(-Effect)
  dat <- flat.dat  %>%
    dplyr::group_by(GeneSymbol) %>%
    dplyr::summarise(MeanLOOCVDiffusionEffect=mean(DiffusionEffect, na.rm=T)) %>%
    ungroup
  dat <- dplyr::left_join(dat,
                          tidyr::spread(flat.dat, loocv, DiffusionEffect),
                          by="GeneSymbol")
  dat
}

#' @noRd
#' @import data.table igraph
#' @importFrom dplyr filter mutate
.knn <- function(hits, loocv.hits, adjm, node.start.count, search.depth, graph)
{
  # TODO diff on loocv.hits
  diffuse.data <- .init.starting.indexes(hits, adjm, node.start.count)
  neighs <- diffusr::nearest.neighbors(
    as.integer(dplyr::filter(diffuse.data$frame, Select==T)$Idx),
    as.matrix(diffuse.data$adjm), search.depth)
  res <- diffuse.data$frame
  neighs
  names(neighs) <- filter(diffuse.data$frame,
                          Idx %in% as.integer(names(neighs)))$GeneSymbol

  neighs <- lapply(neighs, function(e) {
    dplyr::left_join(data.table(Idx=e), res, by="Idx") %>%
      dplyr::select(GeneSymbol, Effect)}
  )
  flat.dat <- do.call(
    "rbind",
    lapply(1:length(neighs),
           function(e) data.table(Start=names(neighs)[e], neighs[[e]])))

  li <- list(data=res, neighbors=flat.dat, graph=graph)
  class(li) <- c("svd.diffused.knn", "svd.diffused")
  return(li)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr left_join mutate select
.init.starting.indexes <- function(hits, adjm, node.start.count)
{
  res  <- dplyr::left_join(data.table::data.table(GeneSymbol=colnames(adjm)),
                           hits, by="GeneSymbol")
  data.table::setDT(res)[is.na(Effect), Effect := 0]
  res <- res %>%
    .[order(-abs(Effect))] %>%
    dplyr::mutate(Idx = 1:.N) %>%
    dplyr::mutate(Select = (Idx <= node.start.count))
  data.table::setDT(res)[is.na(Select), Select := FALSE]
  return(list(frame=res, adjm=adjm))
}
