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
#' @param graph  an weighted adjacency matrix
#' @param r  restart probability of the random if method \code{mrw} is selected
#' @param node.start.count  number of nodes that are used to do the neighbors
#'  search if method \code{neighbors} is selected
#'  If the number of hits in \code{obj} exceeds \code{node.start.count},
#'  then the elements with the highest absolute effects are chosen. If this
#'  behaviour is not desired filter \code{obj} before
#' @param search.depth  how deep should the neighbor search go if method \code{neighbors} is selected
#' @param delete.nodes.on.degree  delete nodes from the graph with a degree of less that \code{delete.nodes.on.degree}
#' @param ...  additional parameters
diffuse <- function(obj, method=c("neighbors", "mrw"), path=NULL, graph=NULL,
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
                                        path=NULL, graph=NULL, r=0.5, node.start.count=25,
                                        search.depth=5, delete.nodes.on.degree=1,
                                        ...)
{
  hits <- obj$gene.hits %>%
    dplyr::select(GeneSymbol, abs(Effect))

  bootstrap.hits <- NULL
  if (any(!is.na(obj$fit$fit$gene.fdrs$FDR)) &
      "MeanBootstrap" %in% colnames(obj$fit$fit$gene.fdrs)) {
    bootstrap.hits <- obj$fit$fit$gene.fdrs %>%
      dplyr::filter(GeneSymbol %in% hits$GeneSymbol) %>%
      dplyr::select(-MeanBootstrap, -Pval, -FDR)
  }

  res   <- .diffuse(hits, bootstrap.hits=bootstrap.hits,
                    path=path, graph=graph,
                    method=match.arg(method),
                    r=r, node.start.count=node.start.count,
                    search.depth=search.depth,
                    delete.nodes.on.degree=delete.nodes.on.degree)
  invisible(res)
}

#' @export
#' @import data.table igraph
#' @importFrom dplyr filter select rename
#' @method diffuse data.table
diffuse.data.table <- function(obj, method=c("neighbors", "mrw"),
                               path=NULL, graph=NULL, r=0.5, node.start.count=25,
                               search.depth=5, delete.nodes.on.degree=1,
                               ...)
{
  if ("Readout" %in% colnames(obj))
  {
    hits <- obj %>%
      dplyr::select(GeneSymbol, Readout) %>%
      dplyr::rename(Effect=Readout)
  }
  else if ("MeanEffect" %in% colnames(obj))
  {
    hits <- obj %>%
      dplyr::select(GeneSymbol, MeanEffect) %>%
      dplyr::rename(Effect=MeanEffect)
  }
  else if ("Effect" %in% colnames(obj))
  {
    hits <- obj %>% dplyr::select(GeneSymbol, Effect)
  }
  else
  {
    stop("Neither 'Readout' nor 'Effect' found in colnames")
  }
  hits <- hits %>% dplyr::mutate(Effect=abs(Effect))
  res   <- .diffuse(hits=hits,
                    bootstrap.hits=NULL,
                    path=path, graph=graph,
                    method=match.arg(method),
                    r=r, node.start.count=node.start.count,
                    search.depth=search.depth,
                    delete.nodes.on.degree=delete.nodes.on.degree)
  res
}

#' @noRd
#' @import data.table igraph
.diffuse <- function(hits, bootstrap.hits, path, graph, method, r,
                     node.start.count, search.depth, delete.nodes.on.degree)
{
  graph <- .read.graph(path=path, graph=graph)
  graph <- igraph::delete.vertices(
    graph, igraph::V(graph)[
      igraph::degree(graph) <= delete.nodes.on.degree  ])
  adjm  <- igraph::get.adjacency(graph, attr="weight")
  l <- switch(method,
         "neighbors"= .knn(hits, bootstrap.hits, adjm, node.start.count, search.depth, graph),
         "mrw"      = .mrw(hits, bootstrap.hits, adjm, r, graph),
         stop("No suitable method found"))
  l
}

#' @noRd
#' @import data.table
#' @importFrom diffusr random.walk
.mrw <- function(hits, bootstrap.hits, adjm, r, graph)
{
  diffuse.data <- .do.mrw(hits, adjm, r)
  if (!is.null(bootstrap.hits))
  {
    pvals <- .significance.mrw(bootstrap.hits, adjm, r)
    res   <- dplyr::left_join(diffuse.data$frame, pvals, by="GeneSymbol")
  }
  else
  {
    res  <- diffuse.data$frame
  }
  li  <- list(diffusion=dplyr::select(res, GeneSymbol, Effect, DiffusionEffect),
              diffusion.model=res,
              lmm.hits=hits,
              graph=graph)
  class(li) <- c("svd.diffused.mrw", "svd.diffused")
  return(li)
}

#' @noRd
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
#' @import foreach
#' @import parallel
#' @import doParallel
.significance.mrw <- function(bootstrap.hits, adjm, r)
{
  boot.g <- tidyr::gather(bootstrap.hits, Boot, Effect, -GeneSymbol) %>%
    as.data.table
  li <- list()
  foreach::foreach(lo=unique(boot.g$Boot)) %do%
  {
    hits <- dplyr::filter(boot.g, Boot==lo)
    hits <- dplyr::select(hits, -Boot)
    dd.lo <- .do.mrw(hits, adjm, r)
    dd.lo$frame$boot <- lo
    li[[lo]] <- dd.lo$frame
  }
  # TODO: how to do hypothesis test here?
  # -> calculate means -> calculate alphas -> calculate variance -> calculate alpha_0
  flat.dat <- do.call("rbind", lapply(li, function(e) e)) %>% dplyr::select(-Effect)
  dat <- flat.dat  %>%
    dplyr::group_by(GeneSymbol) %>%
    dplyr::summarise(MeanBootstrapDiffusionEffect=mean(DiffusionEffect, na.rm=T)) %>%
    ungroup
  dat <- dplyr::left_join(dat,
                          tidyr::spread(flat.dat, boot, DiffusionEffect),
                          by="GeneSymbol")
  dat
}

#' @noRd
#' @import data.table igraph
#' @importFrom dplyr filter mutate
.knn <- function(hits, bootstrap.hits, adjm, node.start.count, search.depth, graph)
{
  # TODO diff on loocv.hits
  diffuse.data <- .init.starting.indexes(hits, adjm, node.start.count)
  neighs <- diffusr::nearest.neighbors(
    as.integer(dplyr::filter(diffuse.data$frame, Select==T)$Idx),
    as.matrix(diffuse.data$adjm), search.depth)
  res <- diffuse.data$frame

  names(neighs) <- (diffuse.data$frame %>%
    dplyr::filter(Idx %in% as.integer(names(neighs))))$GeneSymbol
  # Add the effects of the neighbors
  neighs <- lapply(neighs, function(e) {
    dplyr::left_join(data.table(Idx=e), res, by="Idx") %>%
      dplyr::select(GeneSymbol, Effect)}
  )

  flat.dat <- do.call(
    "rbind",
    lapply(1:length(neighs),
           function(e) data.table(Start=names(neighs)[e], neighs[[e]]))
  )

  genes.found <- dplyr::select(flat.dat, Start, GeneSymbol) %>%
    unlist %>% unname %>% unique
  genes.f.table <- data.table::data.table(GeneSymbol=genes.found) %>%
    dplyr::left_join(res, by="GeneSymbol")

  li <- list(data=res,
             neighbors=flat.dat,
             graph=graph,
             genes.found=genes.f.table)
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
  data.table::setDT(res)[Effect == 0  , Select := FALSE]
  data.table::setDT(res)[is.na(Select), Select := FALSE]
  return(list(frame=res, adjm=adjm))
}
