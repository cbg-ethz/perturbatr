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
#' @importFrom graphics hist
#' @method hist svd.plate
hist.svd.plate <- function(x,  ...)
{
  .hist(x, ...)
}

#' @export
#' @import data.table
#' @importFrom dplyr filter
#' @method hist svd.raw
hist.svd.raw <- function(x,  ...)
{
  ret <- dplyr::filter(x, ReadoutClass=="Readout")
  .hist(ret, ...)
}

#' @export
#' @import data.table
#' @method hist svd.data
hist.svd.data <- function(x,  ...)
{
  .hist(x, ...)
}


#' @noRd
#' @import ggplot2
#' @import data.table
.hist <- function(x,  ...)
{
  df <- data.frame(x=x$Readout)
  pl <- ggplot2::ggplot(df) +
    ggplot2::geom_histogram(aes(x=x, y=..density..),
                            alpha=.15, position="identity", bins=42) +
    ggplot2::geom_density(aes(x)) +
    theme_bw() +
    theme(text = element_text(size = 14, family = "Helvetica")) +
    xlab("Readout") +
    ylab("Density")
  pl
}
