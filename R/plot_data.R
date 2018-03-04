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



#' @title Plot perturbation data
#'
#' @description Creates a barplot of replicate and gene counts of a
#'  \code{PerturbationData} object.
#'
#' @method plot PerturbationData
#' @export
#'
#' @import data.table
#' @import ggplot2
#' @importFrom dplyr summarize
#' @importFrom dplyr group_by
#' @importFrom tidyr gather
#' @importFrom scales pretty_breaks
#'
#' @param x  the object to plot
#' @param size  size of letters
#'
#' @return  returns a plot object
#'
plot.PerturbationData <- function(x, size=10)
{
  dat <- dataSet(x)
  if (dataType(x) == .dataTypes()$RAW && "ReadoutClass" %in% colnames(dat))
    dat <- dplyr::filter(dat, ReadoutClass="Readout")

  pl <-
    dplyr::group_by(dat, Condition) %>%
    dplyr::summarize(Replicates = length(unique(Replicate)),
                     Genes      = length(unique(GeneSymbol))) %>%
    tidyr::gather(Type, Count, Replicates, Genes) %>%
    dplyr::mutate(Count = as.integer(Count)) %>%
    ggplot2::ggplot(ggplot2::aes(x=Condition, y = Count)) +
    ggplot2::geom_bar(ggplot2::aes(fill=Condition), stat="identity") +
    ggplot2::scale_fill_grey(start=0.3) +
    ggplot2::scale_y_continuous(breaks=scales::pretty_breaks(5)) +
    ggplot2::facet_grid(Type ~ ., scales='free_y') +
    ggplot2::geom_text(ggplot2::aes(label = Count, y = Count),
                                    size = floor(size/3), vjust=0) +
    ggplot2::theme_minimal() +
    ggplot2::theme(strip.text      = ggplot2::element_text(size = size),
                   text            = ggplot2::element_text(size = size),
                   panel.spacing.y = ggplot2::unit(2, "lines")) +
    ggplot2::guides(fill=FALSE)

  pl
}
