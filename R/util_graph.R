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
#' @importFrom igraph graph.data.frame
#' @importFrom utils read.csv
.read.graph <- function(pth)
{
  tab <- utils::read.csv(pth, sep="\t", header=T)
  gra <- igraph::graph.data.frame(tab, directed=F)
  if (ncol(tab) == 3)
  {
    if (is.null(E(gra)$weight))
      stop("Third column sould be 'weight'")
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
