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

#' Plot the effects matrix of a PMM analysis
#'
#' @export
#' @import data.table
#'
#' @param x  an svd.prioritized.pmm object
#' @param ...  additional parameters
show.effect.matrix <- function(x, ...)
{
  UseMethod("show.effect.matrix")
}

#' @export
#' @import data.table
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tidyr gather
#' @method show.effect.matrix svd.prioritized.pmm
show.effect.matrix.svd.prioritized.pmm <- function(x, ...)
{
  # TODO: change for other prioritization than abs (see: svd.prioritize.pmm)
  effect.matrices <- effect.matrices(x)
  v <- effect.matrices$gene.effects %>% .[order(GeneSymbol, decreasing=T)]
  v$GeneSymbol <- factor(v$GeneSymbol, levels=v$GeneSymbol)
  v$Pathogen <- "All pathogens"

  LDcolors <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
  pl1 <-
    ggplot2::ggplot(v, aes(Virus, GeneSymbol)) +
    ggplot2::geom_tile(aes(fill=Effect), colour=LDcolors[1]) +
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::scale_y_discrete(expand = c(0,0)) +
    ggplot2::scale_fill_gradient2(low=LDcolors[1], high=LDcolors[11],
                                  na.value=LDcolors[6],
                                  name="Gene effect") +
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size = 14, family = "Helvetica"),
                   aspect.ratio=4,
                   axis.text.x=element_text(angle=45,  hjust = 1, size=10),
                   axis.text.y=element_text(size=10),
                   axis.title=element_blank(), axis.ticks=element_blank())
  v <-  effect.matrices$gene.pathogen.hits %>% gather(GeneSymbol)
  colnames(v) <- c("GeneSymbol", "Virus", "Effect")
  v$GeneSymbol <- factor(v$GeneSymbol, levels=rev(unique(v$GeneSymbol)))
  pl2 <-
    ggplot2::ggplot(v, aes(GeneSymbol, Virus)) +
    ggplot2::geom_tile(aes(fill = Effect), colour=LDcolors[1]) +
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::scale_y_discrete(expand = c(0,0)) +
    ggplot2::scale_fill_gradient2(low=LDcolors[1], high=LDcolors[11],
                                  na.value=LDcolors[6],
                                  name="Gene virus\neffect") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size = 8, family = "Helvetica"),
                   aspect.ratio=2,
                   axis.text.x=element_text(angle=45,  hjust = 1, size=9),
                   axis.text.y=element_text(size=9),
                   axis.title=element_blank(),
                   axis.ticks=element_blank())
  return(list(gene.effects.plot=pl1, gene.pathogen.effects.plot=pl2))
}
