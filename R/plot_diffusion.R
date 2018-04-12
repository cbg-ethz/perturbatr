# perturbatr: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2018 Simon Dirmeier
#
# This file is part of perturbatr
#
# perturbatr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# perturbatr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with perturbatr. If not, see <http://www.gnu.org/licenses/>.


#' @include class_data.R
#' @include class_analysed.R


#' @title Plot a \code{NetworkAnalysedPerturbationData} object
#'
#' @description Creates a table of the gene ranking of a
#' \code{NetworkAnalysedPerturbationData} object
#'
#' @method plot NetworkAnalysedPerturbationData
#' @export
#'
#' @import tibble
#' @import grid
#' @importFrom rlang .data
#'
#' @param x  a \code{NetworkAnalysedPerturbationData} object
#' @param size  size of letters
#' @param main  title of the plot
#' @param ...  additional parameters
#'
#' @return returns a table if the first \code{cnt} highest ranked genes
#'
plot.NetworkAnalysedPerturbationData <- function(x, size=10, main="", ...)
{
  geef <- geneEffects(x)
  geef <- geef[order(-geef$DiffusionEffect), ]
  geef <- geef[1:25,]
  .plot.bars.diff(geef, size, main)
}

#' @noRd
#' @importFrom rlang .data
.plot.bars.diff <- function(x, size, main)
{
  x$GeneSymbol <- base::factor(x$GeneSymbol, levels=rev(unique(x$GeneSymbol)))
  pl <-
    ggplot2::ggplot(x) +
    ggplot2::geom_bar(
      ggplot2::aes(x$GeneSymbol, abs(x$DiffusionEffect)), fill="darkgrey",
      stat="identity") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = size - 2,
                                                       family = "Helvetica"),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   legend.position = "bottom",
                   axis.ticks   = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = size,
                                                family = "Helvetica"),
                   axis.text.x = ggplot2::element_text(
                     size = size - 2,
                     family = "Helvetica")) +
    ggplot2::coord_flip() +
    ggplot2::ggtitle(main)

  pl
}

