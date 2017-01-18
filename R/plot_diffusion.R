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
#' @method plot svd.diffused.mrw
plot.svd.diffused.mrw <- function(x, y, ...)
{
  pars <- list(...)
  sz          <- ifelse(methods::hasArg(size), pars$size, -1)
  show.labels <- ifelse(methods::hasArg(size), pars$size, -1)

  obj <- x$graph.info$graph
  igraph::V(obj)$size = igraph::degree(obj)
  deg <- igraph::degree(obj)
  size <- deg
  size[deg < 3] <- 15
  size[deg >= 3] <- 20
  size[deg > 5] <- 25
  ad <- igraph::get.adjacency(obj)
  ad[ad >= 1] <- 1
  obj <- igraph::graph_from_adjacency_matrix(ad, mode="undirected")
  blue.genes <-
    igraph::V(x$graph.info$graph)[which(igraph::V(x$graph.info$graph)$color == "lightblue")]
  orange.genes <-
    igraph::V(x$graph.info$graph)[which(igraph::V(x$graph.info$graph)$color == "orange")]
  igraph::V(obj)$color[V(obj) %in% blue.genes] <- "lightblue"
  igraph::V(obj)$color[V(obj) %in% orange.genes] <- "orange"
  igraph::E(obj)$width <- 2
  graphics::plot.new()
  op <- par(family = "Helvetica", font=2)
  if (sz != -1) size <- rep(sz, length(size))
  graphics::plot(obj, vertex.size=size,layout =  layout.kamada.kawai,
                 vertex.label.family="Helvetica", vertex.label.font=2,
                 edge.curved=-.01)
  graphics::legend("topright",
                   legend=c("Linear mixed model", "Diffusion"),
                   col=x$graph.info$colors,
                   pch=19, cex=1.05)
  par(op)
}
