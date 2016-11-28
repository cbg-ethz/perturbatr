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
#' @param obj  an svd.data object
#' @param method  method that should be used for diffusion
#' \itemize{
#'   \item{neighbors }{ just looks at the neighbors :)}
#'   \item{mrw }{ do a Markov random walk}
#' }
#' @param path  path to the network
#' @param ...  additional parameters
diffuse <- function(obj, method=c("neighbors", "mrw"), path, ...)
{
  UseMethod("diffuse")
}

#' @noRd
#' @export
#' @import data.table igraph
diffuse.svd.prioritized.pmm <-
function(obj, method=c("neighbors", "mrw"), path, ...)
{
  if (!file.exists(path))
    stop(paste("Can't find: ", path, "!", sep=""))
  hits <- obj$gene.effect.hits %>%
    dplyr::select(GeneSymbol, abs(Effect))
  res  <- .diffuse(hits, path, match.arg(method), ...)
  invisible(res)
}

#' @noRd
#' @import data.table
.diffuse <- function(hits, path, method, ...)
{
  graph <- .read.graph(path)

  f <- switch(method,
         "neighbors"= .knn, "mrw" = .mrw,
         stop("No suitable method found"))
  f(hits, graph, ...)
}

#' @noRd
#' @import data.table
#' @importFrom diffusr random.walk
.mrw <- function(hits, graph, r=0.5, ...)
{
  diffuse.data <- .init.starting.distribution(hits, graph)
  mrw <- diffusr::random.walk(abs(diffuse.data$frame$Effect),
                              diffuse.data$adjm, .5)
  diffuse.data$frame$DiffusionEffect <- mrw
  res <- diffuse.data$frame
  class(res) <- c("svd.diffused.mrw", "svd.diffused", class(res))
  return(res)
}


#' @noRd
#' @import data.table igraph
#' @import foreach parallel doParallel
#' @importFrom iterators iter
#' @importFrom dplyr filter
#' @importFrom methods hasArg
.knn <- function(hits, adjm, ...)
{
  pars <- list(...)
  c <- ifelse(hasArg(count), pars$count, 1)
  phs <- obj$gene.pathogen.effect.hits
  vir <- unique(phs$Virus)
  all.edges  <- igraph::get.edgelist(graph)
  #cl <- parallel::makeCluster(parallel::detectCores() - 2)
  #doParallel::registerDoParallel(cl)
  neighbors <- foreach::foreach(v=iterators::iter(vir), .combine=rbind) %do%
  {
      virgen <- dplyr::filter(phs, Virus==v)$GeneSymbol
      idxs   <- which(all.edges[,1] %in% virgen | all.edges[,2] %in% virgen)
      fr <- t(apply(all.edges[idxs, ], 1, function(e) sort(e)))
      chosen.edges <- data.table(Virus=v,
                                    Gene1=c(fr[, 1]),
                                    Gene2=c(fr[, 2])) %>%
        unique
      chosen.edges
  }
 # parallel::stopCluster(cl)
  rel <- dplyr::group_by(neighbors, Gene1, Gene2) %>%
    dplyr::summarise(Count=n()) %>%
    dplyr::filter(Count >= c) %>%
    as.data.frame
  li <- unique(c(rel$Gene1, rel$Gene2))
  res <- dplyr::filter(neighbors, Gene1 %in% li | Gene2 %in% li ) %>%
    as.data.frame
  nodes <- data.table(Node=unique(c(res[, 2], res[, 3]))) %>%
    dplyr::mutate(Color=ifelse(Node %in% phs$GeneSymbol,
                               "lightblue", "orange")) %>%
    dplyr::mutate(FromLMM=ifelse(Node %in% phs$GeneSymbol, 1, 0)) %>%
    as.data.frame
  edges <- res[, c(2, 3)]
  res.gr <- igraph::graph.data.frame(edges, directed=F, vertices=nodes)
  igraph::V(res.gr)$color <- igraph::V(res.gr)$Color
  graph.info <- list(graph=res.gr,
                     colors=c("lightblue", "orange"),
                     legend=c("Linear mixed model", "Diffusion"),
                     type="1-NN",
                     tresh=2)
  list(hits=res, graph.info=graph.info)
}

#' @noRd
#' @import data.table
#' @importFrom igraph get.adjacency
#' @importFrom dplyr left_join
.init.starting.distribution <- function(hits, graph)
{
  adjm <-  igraph::get.adjacency(graph, attr="weight")
  res  <- dplyr::left_join(data.table::data.table(GeneSymbol=colnames(adjm)),
                           hits)
  data.table::setDT(res)[is.na(Effect), Effect := 0]
  return(list(frame=res, adjm=adjm))
}
