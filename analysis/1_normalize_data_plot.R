library(data.table)
library(dtplyr)
library(dplyr)
library(tidyr)
library(lme4)
library(optparse)
library(knockout)
library(ggplot2)
library(hashmap)
library(ggthemr)
library(hrbrthemes)
library(viridis)
library(cowplot)

ggthemr("fresh", "scientific")

plot.data <- function(x)
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
    ggplot2::scale_y_continuous(breaks=scales::pretty_breaks(5)) +
    ggplot2::facet_grid(Type ~ Screen, scales='free_y') +
    ggplot2::geom_text(ggplot2::aes(label = Count, y = Count), size = floor(20/3), vjust=0) +
    ggplot2::theme_bw() +
    hrbrthemes::theme_ipsum_rc(base_family="Helvetica") +
    ggplot2::theme(text            = ggplot2::element_text(size = 20),
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
  qual <- x@.data
  qual <- dplyr::group_by(qual, Virus, Screen, Library,
                          ScreenType, ReadoutType,
                          Replicate, Plate) %>%
    dplyr::mutate(Plate = .GRP) %>%
    ungroup
  df <- qual %>%
    dplyr::select(Virus, Screen, Readout, Plate) %>%
    dplyr::mutate(Plate=as.factor(Plate)) %>%
    dplyr::filter(!is.na(Readout), !is.infinite(Readout))


  ggplot2::ggplot(df, ggplot2::aes(Plate, Readout)) +
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

  path <- "./"
  out.dir <- "./plots"

  raw.data.file  <- paste(path, "data/rnai_screen_raw.rds", sep="/")
  norm.data.file <- paste(path, "data/rnai_screen_normalized.rds", sep="/")

  raw.dat  <- readRDS(raw.data.file)
  norm.dat <- readRDS(norm.data.file)

  pl <- plot.data(norm.dat)


  ggsave(
    filename = paste0(out.dir, "/",  "data_overview-", ".eps"),
    plot = pl,
    width = 10,
    height = 10)

  pl.raw  <- knockout::filter(raw.dat, Virus=="HCV", Screen=="Kinome") %>%
    plot.quality
  pl.norm  <- knockout::filter(norm.dat, Virus=="HCV", Screen=="Kinome") %>%
    plot.quality + ylab("")

  pl.gr <- plot_grid(pl.raw, pl.norm, align="h", labels=c("(a)", "(b)"))

  ggsave(
    filename = paste0(out.dir, "/", "comparison_plate_readouts", ."eps")
    plot = pl.gr,
    width = 12,
    height = 6
  )

}

run()
