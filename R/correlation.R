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

#' Calculate the correlation between two data-sets
#'
#' @export
#' @import data.table
#'
#' @param x  a random object
#' @param y  a random object of the same type (optional)
#' @param method  a character string indicating which correlation coefficient
#' (\emph{pearson}, \emph{kendall} or \emph{spearman})
#' @param ...  additional parameters
correlation <- function(x, y, method = c("pearson", "kendall", "spearman"), ...)
{
  UseMethod("correlation")
}

#' @export
#' @import data.table
#' @importFrom dplyr filter
#' @importFrom stats cor
#' @method correlation svd.replicates
correlation.svd.replicates <- function(x, y,
                                       method = c("pearson", "kendall", "spearman"),
                                       ...)
{
  method <- match.arg(method)
  len <- length(x)
  gri <- expand.grid(1:len, 1:len)
  ret <- do.call(
    "rbind",
    apply(gri, 1, function (gr)
    {
      val1 <- unname(gr[1])
      val2 <- unname(gr[2])
      x1 <- unname(x[[val1]]$Readout)
      x2 <- unname(x[[val2]]$Readout)
      idx <- which(!is.na(x1) & !is.na(x2))
      df <- data.table::data.table(
        Var1=val1, Var2=val2,
        Cor=stats::cor(x1[idx], x2[idx], method=method))
      df
    })
  )
  ret <- data.table::dcast(ret, Var1 ~ Var2, value.var="Cor") %>%
    dplyr::select(-Var1) %>%
    as.matrix(ret)
  colnames(ret) <- rownames(ret) <- paste("Replicate", 1:len, sep="_")
  ret
}
