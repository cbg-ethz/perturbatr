#!/usr/bin/env Rscript

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


library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(perturbatr)


plot.data <- function(x)
{
  numb.frame <-
    dplyr::group_by(x, Condition, Screen) %>%
    dplyr::summarize(Replicates = length(unique(Replicate)),
                     Genes      = length(unique(GeneSymbol))) %>%
    tidyr::gather(Type, Count, Replicates, Genes)
  numb.frame$Count <- as.integer(numb.frame$Count)

  pl <-
    ggplot2::ggplot(numb.frame, ggplot2::aes(x=Condition, y = Count)) +
    ggplot2::geom_bar(ggplot2::aes(fill=Condition), stat="identity") +
    ggplot2::scale_y_continuous(breaks=scales::pretty_breaks(5)) +
    ggplot2::facet_grid(Type ~ Screen, scales='free_y') +
    ggplot2::geom_text(ggplot2::aes(label = Count, y = Count), size = floor(20/3), vjust=0) +
    ggplot2::theme_bw() +
    hrbrthemes::theme_ipsum_rc(base_family="Helvetica") +
    ggplot2::theme(text = ggplot2::element_text(size = 20),
                   axis.text.x=element_blank(),
                   axis.text.y=element_text(size=14),
                   axis.title.x=element_blank(),
                   axis.title.y=element_text(size=16),
                   legend.text=element_text(size=14),
                   plot.title = element_text(hjust = 0.5),
                   strip.text.x = element_text(size = 14, hjust=.1),
                   strip.text.y = element_text(size = 14),
                   panel.spacing.y = ggplot2::unit(2, "lines"))
}


plot.quality <- function(x)
{
  qual <- x
  qual <- dplyr::group_by(qual, Condition, Screen, Library,
                          ScreenType, ReadoutType,
                          Replicate, Plate) %>%
    {dplyr::mutate(ungroup(.), grp = group_indices(.)) }%>%
    ungroup()
  df <- qual %>%
    dplyr::select(Condition, Screen, Readout, grp) %>%
    dplyr::mutate(grp=as.factor(grp)) %>%
    dplyr::filter(!is.na(Readout), !is.infinite(Readout))


  ggplot2::ggplot(df, ggplot2::aes(grp, Readout)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::theme(aspect.ratio=.75,
      text = ggplot2::element_text(size = 20),
                   axis.text.x=element_blank(),
                   axis.text.y=element_text(size=14),
                   axis.title.x=element_text(size=16),
                   axis.title.y=element_text(size=16),
                   legend.text=element_text(size=14)) +
    xlab("Plate")

}


run <- function()
{

  data.dir <- "./data"
  out.dir  <- "./plots"

  raw.data.file  <- paste(data.dir, "rnai_screen_raw.rds", sep="/")
  norm.data.file <- paste(data.dir, "rnai_screen_normalized_2.rds", sep="/")

  raw.dat  <- readRDS(raw.data.file)
  norm.dat <- readRDS(norm.data.file)

  pl <- plot.data(norm.dat)

  ggsave(
    filename = paste0(out.dir, "/",  "data_overview", ".eps"),
    plot = pl,
    width = 10,
    height = 10)

  pl.raw  <- dplyr::filter(raw.dat, Condition == "HCV", Screen=="Kinome") %>%
    plot.quality()
  pl.norm  <- dplyr::filter(norm.dat, Condition == "HCV", Screen=="Kinome") %>%
    plot.quality + ylab("")
  pl.gr <- cowplot::plot_grid(pl.raw, pl.norm, align="h", labels=c("(a)", "(b)"))

  ggsave(
    filename = paste0(out.dir, "/", "comparison_plate_readouts", ".eps"),
    plot = pl.gr,
    width = 12,
    height = 6)
}

run()
