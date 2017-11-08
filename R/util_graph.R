# knockdown: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2015 - 2016 Simon Dirmeier
#
# This file is part of knockdown
#
# knockdown is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# knockdown is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with knockdown. If not, see <http://www.gnu.org/licenses/>.


#' @noRd
#' @importFrom utils read.csv
#' @importFrom igraph graph_from_adjacency_matrix graph.data.frame
.read.graph <- function(path, graph)
{
  if (all(is.null(c(path, graph))))
  {
    stop("Please provide either a graph or the file-to a graph!")
  }
  if (!is.null(path))
  {
    if (!file.exists(path)) stop("'path' does not exist!")
    tab <- utils::read.csv(path, sep="\t", header=TRUE)
    gra <- igraph::graph.data.frame(tab, directed=FALSE)
    if (ncol(tab) == 3)
    {
      if (is.null(igraph::E(gra)$weight))
        stop("Third column sould be 'weight'")
    }
  }
  else if(!is.null(graph))
  {
    gra <- igraph::graph_from_adjacency_matrix(
      graph, weighted=TRUE, mode="undirected")
  }
  else
  {
    stop("Something went wrong with graph-reading.")
  }
  gra
}
