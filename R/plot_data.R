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

#' Plot a knockout dataset
#'
#' @export
#' @import data.table
#' @method plot knockout.raw.data
#' @param x  the object to plot
#' @param ...  additional parameters
plot.knockout.raw.data <- function(x, size=10, ...)
{
  x@.data <- dplyr::filter(x@.data, ReadoutClass=="Readout")
  plot.knockout.normalized.data(x, size, ...)
}

#' Plot a knockout data-set
#'
#' @export
#' @method plot knockout.normalized.data
#' @import ggplot2
#' @import data.table
#' @importFrom dplyr summarize
#' @importFrom dplyr group_by
#' @importFrom tidyr gather
#' @param x  the object to plot
#' @param ...  additional parameters
plot.knockout.normalized.data <- function(x, size, ...)
{
  numb.frame <-
    dplyr::group_by(x@.data, Virus, Screen) %>%
    dplyr::summarize(Replicates = length(unique(Replicate)),
                     Genes      = length(unique(GeneSymbol))) %>%
    tidyr::gather(Type, Count, Replicates, Genes)

  numb.frame$Count <- as.integer(numb.frame$Count)
  pl <-
    ggplot2::ggplot(numb.frame, ggplot2::aes(x=Virus, y = Count)) +
    ggplot2::geom_bar(ggplot2::aes(fill=Virus), stat="identity") +
    ggplot2::facet_grid(Type ~ Screen, scales='free_y') +
    ggplot2::scale_fill_brewer(palette="Spectral") +
    ggplot2::geom_text(ggplot2::aes(label = Count, y = Count), size = floor(size/3), vjust=0) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = size),
                   text       = ggplot2::element_text(size = size))

  pl
}
