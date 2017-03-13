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

#' @noRd
#' @import igraph
#' @importFrom utils read.csv
.read.graph <- function(path, graph)
{
  if (all(is.null(c(path, graph)))) {
    stop("Please provide either a graph or the file-to a graph!")
  }

  if (!is.null(path))
  {
    if (!file.exists(path)) stop("'path' does not exist!")
    tab <- utils::read.csv(path, sep="\t", header=T)
    gra <- igraph::graph.data.frame(tab, directed=F)
    if (ncol(tab) == 3)
    {
      if (is.null(igraph::E(gra)$weight))
        stop("Third column sould be 'weight'")
    }
  }
  else if(!is.null(graph))
  {
    gra <- graph_from_adjacency_matrix(graph, weighted=T, mode="undirected")
  }
  else
  {
    stop("Something went wrong with graph-reading.")
  }
  gra
}

#' @noRd
#' @import Matrix
#' @importFrom assertthat assert_that
.stoch.col.norm <- function(m)
{
  col.sums  <- Matrix::colSums(m)
  zero.cols <- Matrix::which(col.sums == 0)
  if (length(zero.cols) != 0)
  {
    m[,zero.cols]       <- 1
    col.sums[zero.cols] <- nrow(m)
  }
  ret <- sweep(m, 2L, col.sums, "/", check.margin = FALSE)
  ret.col.sums <- unname(Matrix::colSums(ret))
  assertthat::assert_that(all(ret.col.sums <= 1.001))
  assertthat::assert_that(all(ret.col.sums >= 0.999))
  assertthat::assert_that(all(ret >= 0) & all(!is.nan(ret@x)))
  invisible(ret)
}

#' Get the neighbors of a node including all edges between them
#'
#' @param gene  the gene for which the neighborhood is searched
#' @param graph  the graph
#' @import igraph
.induced.subgraph <- function(gene, graph)
{
  stop("not yet imlemented")
  fr  <- t(apply(edge.list, 1, function(e) sort(e)))

  edge.list <- igraph::get.edgelist(graph)
  idxs <- edge.list[,1] == gene | edge.list[,2] == gene
  genes <- edge.list[idxs, ]
  genes
}
