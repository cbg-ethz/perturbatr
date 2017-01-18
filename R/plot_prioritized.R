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
#' @import ggplot2
#' @import data.table
#' @importFrom dplyr filter
#' @method plot svd.prioritized.pmm
plot.svd.prioritized.pmm <- function(x, y, ...)
{
  pl <- .plot.svd.prioritized.pmm(x$gene.hits, main="Gene effects", ...)
  pl2 <-
    .plot.svd.prioritized.pmm(x$gene.pathogen.hits, main="Gene-virus effects", ...) +
    ggplot2::facet_wrap( ~ Virus, ncol=ceiling(length(unique(gen.pat$Virus))/2))
  return(list(pl, pl2))
}

#' @noRd
#' @import data.table
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom methods hasArg
.plot.svd.prioritized.pmm <- function(x, main, ...)
{
  pars <- list(...)
  size <- ifelse(methods::hasArg(size), pars$size, 10)
  if ("Virus" %in% colnames(x))
  {
    x <- filter(x, Control == 0) %>%
      .[order(abs(Effect), decreasing=T), .SD[1:25], by=Virus] %>%
      dplyr::filter(!is.na(GeneSymbol))
  }
  else
  {
    x <- x[order(abs(Effect), decreasing=T), .SD[1:25]] %>%
      dplyr::filter(!is.na(GeneSymbol), !is.na(Effect))
  }
  x.pos.range <- max(abs(x$Effect))
  x.lim  <- c(-x.pos.range, x.pos.range) + c(-x.pos.range, x.pos.range)/5
  pl <-
    ggplot2::ggplot(x) +
    ggplot2::geom_bar(aes(x=GeneSymbol,y=abs(Effect), fill=Effect),
                      stat="identity") +
    ggplot2::scale_fill_distiller(palette="Spectral", limits=x.lim) +
    ylab("Effect") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y=element_blank(),
                   axis.ticks=element_blank(),
                   text = element_text(size = size, family = "Helvetica")) +
    ggplot2::coord_polar() +
    ggplot2::ggtitle(main)

  pl
}

#' @export
#' @import data.table
#' @method plot svd.prioritized.hyper
plot.svd.prioritized.hyper <- function(x, y, ...)
{
  .plot.svd.prioritized(x, ...)
}

#' @export
#' @import data.table
#' @method plot svd.prioritized.tt
plot.svd.prioritized.tt <- function(x, y, ...)
{
  .plot.svd.prioritized(x, ...)
}

#' @import data.table
#' @import ggplot2
#' @importFrom dplyr group_by summarize mutate filter
.plot.svd.prioritized <- function(x, ...)
{
  x <- x[order(abs(MeanEffect), decreasing=T), .SD[1:25]] %>%
    dplyr::filter(!is.na(GeneSymbol), !is.na(MeanEffect))
  x.pos.range <- max(abs(x$MeanEffect))
  x.lim  <- c(-x.pos.range, x.pos.range) + c(-x.pos.range, x.pos.range)/5
  pl <-
    ggplot2::ggplot(x) +
    ggplot2::geom_bar(aes(x=GeneSymbol,y=abs(MeanEffect), fill=MeanEffect),
                      stat="identity") +
    ggplot2::scale_fill_distiller(palette="Spectral", limits=x.lim) +
    ylab("Effect") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y=element_blank(),
                   axis.ticks=element_blank()) +
    ggplot2::coord_polar()
  pl
}
