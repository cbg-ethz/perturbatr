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
#' @import data.table
#' @importFrom RankAggreg RankAggreg
#' @param obj  a n x k matrix where n is the number of lists and k the number of elements in each list
#' @param labels  the names of the k elements
#' @param k  how many elements should be considered
.rankaggreg <- function(obj, labels, k)
{
  ranks <- matrix(NA, nrow(obj), ncol(obj))
  ranks <- t(apply(obj, 2, function(el)
  {
      tab <- data.table::data.table(lab=labels, stat=el) %>%
        .[order(el, decreasing=T)]
      tab$lab
  }))
  ret <- RankAggreg::RankAggreg(ranks[,1:100],
                                method="CE",distance="Kendall", k=k)
  ret
}
