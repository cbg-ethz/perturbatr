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
#' @method plot svd.raw
plot.svd.raw <- function(x, y, ...)
{
  x <- dplyr::filter(x, ReadoutClass=="Readout")
  plot.svd.data(x, ...)
}

#' @export
#' @import ggplot2
#' @import data.table
#' @importFrom dplyr summarize
#' @importFrom dplyr group_by
#' @importFrom tidyr gather
#' @method plot svd.data
plot.svd.data <- function(x, y, ...)
{
  numb.frame <-
    dplyr::group_by(x, Virus, Screen) %>%
    dplyr::summarize(Replicates=length(unique(Replicate)),
                     Genes=length(unique(GeneSymbol))) %>%
    tidyr::gather(Type, Count, Replicates, Genes)
  numb.frame$Count <- as.integer(numb.frame$Count)
  pl <-
    ggplot2::ggplot(numb.frame, aes(x=Virus, y = Count)) +
    ggplot2::geom_bar(aes(fill=Virus), stat="identity") +
    ggplot2::facet_grid(Type ~ Screen, scales='free_y') +
    ggplot2::scale_fill_brewer(palette="Spectral") +
    ggplot2::geom_text(aes(label = Count, y = Count), size = 4, vjust=.25) +
    ggplot2::theme_bw()

  pl
}
