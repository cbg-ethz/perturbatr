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

#' Plot the hit counts of a PMM analysis for every virus as barplots
#'
#' @export
#' @import data.table
#'
#' @param x  an svd.prioritized.pmm object
#' @param ...  additional parameters
show.pathogen.hit.counts <- function(x, ...)
{
  UseMethod("show.pathogen.hit.counts")
}

#' @export
#' @importFrom dplyr group_by summarize mutate
#' @method show.pathogen.hit.counts svd.prioritized.pmm
show.pathogen.hit.counts.svd.prioritized.pmm <- function(x, ...)
{
  obj <- x$gene.pathogen.effect.hits
  single.res <-
    dplyr::group_by(obj, Virus, Sign=sign(Effect)) %>%
    dplyr::summarize(cnt=n()) %>%
    ungroup %>%
    dplyr::mutate(Count=cnt*Sign)

  pl <-
    ggplot2::ggplot(single.res, aes(x=Virus, y=Count, fill=Sign)) +
    ggplot2::geom_bar(stat="identity", position="dodge") +
    ggplot2::scale_fill_distiller(palette="Spectral") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="none",
                   text = element_text(size=22),
                   axis.text.y=element_blank()) +
    ggplot2::geom_hline(yintercept=0) +
    ggplot2::ylab("Count hits") +
    ggplot2::geom_text(aes(x=Virus, y=ifelse(Count>0, Count+1, Count-1),
                           label=abs(Count)), size=5, colour="black")
  pl
}
