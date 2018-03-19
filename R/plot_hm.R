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
#' @include util_enums.R


#' @noRd
#' @export
#' @method plot HMAnalysedPerturbationData
#' @import ggplot2
#' @import tibble
#' @importFrom dplyr filter
#'
#' @param x  the object to be plotted
#' @param size  size of letters
#' @param ...  additional parameters
#'
#' @return returns a list of plots
#'
plot.HMAnalysedPerturbationData <- function(x, size=10, main="", ...)
{
  pl <- try({
    .plot.perturbation.hm.analysed(
      geneHits(x), main=main, size=size)
  })
  pl2 <- try({
    .plot.perturbation.hm.analysed(nestedGeneHits(x), main="", size=size) +
      ggplot2::facet_wrap( ~ Condition,
        ncol=ceiling(length(unique(nestedGeneHits(x)$Condition))/2))
  })
  pl3 <- try({
    .plot.effect.matrices.perturbation.analysed.hm(x, size)
  })

  return(list(gene.effect.barplot        = pl,
              nested.gene.effect.barplot = pl2,
              nested.gene.effect.matrix  = pl3))
}


#' @noRd
#' @import tibble
#' @import ggplot2
#' @importFrom dplyr group_by filter row_number
#' @importFrom methods hasArg
#' @importFrom rlang .data
.plot.perturbation.hm.analysed  <- function(x, main, size, ...)
{
  pars <- list(...)
  x <- dplyr::filter(x, .data$Control == 0)
  x <- x[base::order(-abs(x$Effect)), ]
  if ("Condition" %in% colnames(x))
  {
    x <- dplyr::group_by(x, .data$Condition)
  }
  x <- dplyr::filter(x, row_number() <= 25)
  x$GeneSymbol <- base::factor(x$GeneSymbol, levels=rev(unique(x$GeneSymbol)))
  pl <- .plot.bars(x, size, main, ...)
  pl
}


#' @noRd
#' @import tibble
#' @import ggplot2
#' @importFrom tidyr gather
#' @importFrom rlang .data
.plot.effect.matrices.perturbation.analysed.hm <- function(x, size, ...)
{

  effect.matrices <- effect.matrices(x)
  ge <- effect.matrices$gene.effects
  ge <- dplyr::arrange(ge, desc(abs(.data$Effect)))
  ge <- ge[seq(25), ]

  gpe <- effect.matrices$nested.gene.effects
  gpe <- dplyr::filter(gpe, .data$GeneSymbol %in% ge$GeneSymbol)
  gpe <- tidyr::gather(gpe, "Condition", "Effect", -.data$GeneSymbol)
  gpe$GeneSymbol <- factor(gpe$GeneSymbol, levels=rev(unique(gpe$GeneSymbol)))

  pl <-
    ggplot2::ggplot(gpe, ggplot2::aes(gpe$GeneSymbol, gpe$Condition)) +
    ggplot2::geom_tile(ggplot2::aes(fill = gpe$Effect), colour="black") +
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::scale_y_discrete(expand = c(0,0)) +
    ggplot2::scale_fill_gradient2(low      = colors()$blue,
                                  high     = colors()$red,
                                  na.value = "white",
                                  name     = "Nested gene effect") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::theme(text=ggplot2::element_text(size = size,
                                              family = "Helvetica"),
                   aspect.ratio = 2,
                   axis.text.x  = ggplot2::element_text(angle=45,
                                                        hjust = 1,
                                                        size=size),
                   axis.text.y  = ggplot2::element_text(size=size),
                   axis.title   = ggplot2::element_blank(),
                   axis.ticks   = ggplot2::element_blank())
  pl
}


#' @noRd
#' @importFrom rlang .data
.plot.bars <- function(x, size, main, ...)
{
    ggplot2::ggplot(x) +
    ggplot2::geom_bar(
      ggplot2::aes(x$GeneSymbol, abs(x$Effect), fill=factor(sign(x$Effect))),
      stat="identity") +
    ggplot2::scale_fill_manual("Trend",
      values = c(colors()$red, "grey", colors()$blue),
      limits = c("1",  "0", "-1"),
      labels=c("Positive", "None", "Negative")) +
    ggplot2::scale_x_discrete(labels = rev(sort(x$GeneSymbol)),
                              limits = rev(sort(x$GeneSymbol))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = size - 2,
                                                       family = "Helvetica"),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.ticks=ggplot2::element_blank(),
                   text = ggplot2::element_text(size = size,
                                                family = "Helvetica"),
                   axis.text.x = ggplot2::element_text(
                                              size = size - 2,
                                              family = "Helvetica")) +
    ggplot2::coord_flip() +
    ggplot2::ggtitle(main)
}
