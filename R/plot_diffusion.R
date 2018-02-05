# perturbR: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2018 Simon Dirmeier
#
# This file is part of perturbR
#
# perturbR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# perturbR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with perturbR. If not, see <http://www.gnu.org/licenses/>.


#' @include class_data.R
#' @include class_analysed.R


#' Plot a \code{perturbation.diffusion.analysed} object
#'
#' @method plot perturbation.diffusion.analysed
#'
#' @export
#'
#' @import data.table
#' @import igraph
#' @importFrom graphics plot legend
#' @importFrom methods hasArg
#' @importFrom assertthat assert_that
#'
#' @param x  a \code{perturbation.diffusion.analysed} object
#' @param graph.size  approximate numer of nodes
#' @param show.node.labels  \code{logical} if gene names are shown or note
#' @param ...  addutional params
#'
#' @return returns a plot object
plot.perturbation.diffusion.analysed <- function(x,
                                             graph.size = 20,
                                             show.node.labels = FALSE,
                                             ...)
{
  pars        <- list(...)
  sz          <- ifelse(methods::hasArg(size), pars$size, -1)
  obj         <- x@.graph
  stopifnot(igraph::is.igraph(obj))
  # get the best 100 hits according to their ranking
  v.cnt <- 100
  repeat
  {
    best.100 <- x@.data %>%
      .[order(-DiffusionEffect)] %>%
      .[1:v.cnt] %>%
      dplyr::filter(!is.na(DiffusionEffect))
    # index of edges that are gonna be plotted
    edge.list    <- igraph::get.edgelist(x@.graph)
    idxs         <- .edge.indexes(edge.list, best.100)
    obj          <- igraph::graph.data.frame(.edge.subset(edge.list, idxs),
                                             directed=FALSE)

    if (length(igraph::V(obj)) >= 100) break
    v.cnt <- v.cnt + 1
  }
  # get the connected components
  comps              <- igraph::components(obj)
  # get the genes that are not in the largest component
  non.max.comp.genes <- names(
    which(comps$membership != which.max(comps$csize))
  )
  # remove the genes that are not in the largest component
  obj <- igraph::delete.vertices(obj, non.max.comp.genes)
  # node colors (LMM identified genes are blue, rest orange)
  blue.genes   <- best.100$GeneSymbol[best.100$GeneSymbol %in%
                                        x@.initial.model@.gene.hits$GeneSymbol]
  # set some vis params

  size <- .size(obj)
  igraph::V(obj)$color[size <= 10] <- "#fe9929"
  igraph::V(obj)$color[size == 3]  <- "#fed98e"
  igraph::V(obj)[igraph::V(obj)$name %in% blue.genes] $color <- "blue"
  # igraph::V(obj)$color[igraph::V(obj)$name %in% blue.genes] <- "lightblue"
  igraph::E(obj)$width <- 2

  # plot all the the nodes
  .plot.graph(obj, sz, size, show.node.labels)
}

#' @noRd
.edge.indexes <- function(edge.list, best.100)
{
  v1 <- (edge.list[,1] %in% best.100$GeneSymbol &
           edge.list[,2] %in% best.100$GeneSymbol)
  v2 <- (edge.list[,2] %in% best.100$GeneSymbol &
           edge.list[,1] %in% best.100$GeneSymbol)
  which(v1 | v2)
}

#' @noRd
.size <- function(obj)
{
  deg                 <- igraph::degree(obj)
  size                <- deg
  size[deg <  3]      <- 3
  size[deg >= 3]      <- 5
  size[deg >  5]      <- 10
  unname(size)
}

#' @noRd
#' @import data.table
.edge.subset <- function(edge.list, idxs)
{
  fr  <- t(apply(edge.list[idxs, ], 1, function(e) sort(e)))
  rel <- data.table::data.table(Gene1=c(fr[, 1]),
                                Gene2=c(fr[, 2])) %>%
    unique
  as.data.frame(rel)
}

#' @noRd
#' @importFrom graphics plot par legend
#' @importFrom igraph V layout.kamada.kawai
.plot.graph <- function(obj, sz, size, show.node.labels)
{
  op                  <- graphics::par(family = "Helvetica", font=1)
  if (sz != -1) size  <- rep(sz, length(size))
  if (!show.node.labels)
  {
    igraph::V(obj)$name <- rep(NA, length(igraph::V(obj)$name))
  }
  else
  {
    lbgens <- length(igraph::V(obj)$name[igraph::V(obj)$color == "blue"])
    igraph::V(obj)$name[igraph::V(obj)$color == "blue"] <- rep(NA, lbgens)
  }

  graphics::plot(
    obj,
    vertex.size = size,
    layout = igraph::layout.kamada.kawai,
    edge.curved = -.05
  )
  graphics::legend(
    "topleft",
    box.lty=NULL,
    bty = "n",
    border = FALSE,
    fill = c("blue", "#e34a33", "#fdbb84"),
    legend = c("Host factors identified by LMM",
               "Host factors identified by MRW")
  )
  graphics::par(op)
}
