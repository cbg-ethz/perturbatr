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

#' @export
#' @import data.table
#' @import igraph
#' @importFrom graphics plot legend
#' @importFrom methods hasArg
#' @importFrom assertthat assert_that
#' @method plot svd.diffused.mrw
plot.svd.diffused.mrw <- function(x, y, graph.size=20, ...)
{
  pars        <- list(...)
  sz          <- ifelse(methods::hasArg(size), pars$size, -1)
  obj         <- x$graph
  stopifnot(igraph::is.igraph(obj))
  # get the best 100 hits according to their ranking
  v.cnt <- 100
  repeat
  {
    best.100 <- x$diffusion %>%
      .[order(-DiffusionEffect)] %>%
      .[1:v.cnt] %>%
      dplyr::filter(!is.na(DiffusionEffect))
    # index of edges that are gonna be plotted
    edge.list    <- igraph::get.edgelist(x$graph)
    idxs         <- .edge.indexes(edge.list, best.100)
    obj <- igraph::graph.data.frame(.edge.subset(edge.list, idxs), directed=F)
    if (length(igraph::V(obj)) >= 100) break
    v.cnt <- v.cnt + 1
  }
  # node colors (LMM identified genes are blue, rest orange)
  blue.genes   <- best.100$GeneSymbol[best.100$GeneSymbol %in%
                                        x$lmm.hits$GeneSymbol]
  # set some vis params

  size <- .size(obj)
  igraph::V(obj)$color <- "orange"
  igraph::V(obj)[igraph::V(obj)$name %in% blue.genes] $color <- "blue"
  # igraph::V(obj)$color[igraph::V(obj)$name %in% blue.genes] <- "lightblue"
  igraph::E(obj)$width <- 2
  # plot all the the nodes
  .plot.graph(obj, sz, size)
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
  size[deg <  3]      <- 1
  size[deg >= 3]      <- 3
  size[deg >  5]      <- 5
  size
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
.plot.graph <- function(obj, sz, size)
{
  op                 <- par(family = "Helvetica", font=2)
  if (sz != -1) size <- rep(sz, length(size))
  graphics::plot(obj, vertex.size=size,layout= igraph::layout.kamada.kawai,
                 vertex.label.family="Courier", vertex.label.font=1,
                 vertex.label.cex=.95, edge.curved=-.15)
  par(op)
}
