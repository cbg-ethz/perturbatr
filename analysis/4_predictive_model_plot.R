#!/usr/bin/env Rscript

library(dplyr)
library(data.table)
library(tidyr)
library(lme4)
library(optparse)
library(knockdown)
library(ggplot2)
library(ggthemr)
library(hrbrthemes)
library(viridis)
library(cowplot)

ggthemr("fresh", "scientific")

jaccard <- function(s1, s2) { length(intersect(s1, s2)) / length(union(s1, s2)) }

bio.pred <- function(fls.bio)
{

  stab <- readRDS(fls.bio)
  stab <- stab[-1]

  df <- rbindlist(lapply(stab, function(el) {
    data.table(Virus=length(unlist(stringr::str_split(el$Virus, "_"))),
               el$sse)}))

  df <- df %>% dplyr::filter(Sampling=="CV")

  df$Model[df$Model == "LMM"] <- "knockdown"
  df$Virus <- paste(df$Virus, "viruses")

  p <-
    ggplot2::ggplot(df) +
    geom_boxplot(aes(x     = Model,
                     y     = MSE,
                     fill  = Model), width=.5) +
    ggplot2::facet_grid(Virus ~  .) +
    ggplot2::scale_fill_manual(values =ggthemr:::palettes$fresh$swatch[c(2,4)]) +
    ggplot2::theme_bw() +
    hrbrthemes::theme_ipsum_rc(base_family="Helvetica") +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_text(size=14),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=16),
          legend.text=element_text(size=14),
          text = ggplot2::element_text(size = 20),
          strip.text.y = element_text(size = 14)) +
    ylab("Mean residual sum of squares")

  p
}

syn.pred <- function(fls.stn)
{
    df <- rbindlist(lapply(fls.stn, function (f) {
      fld <- readRDS(f)
      me <- stringr::str_match(f, ".*genemean_(.+).rds")[2] %>% as.numeric
      dff <- NULL
      ix <- grep("bt$", names(fld$benchmark))
      for (i in ix)
      {
        curd <- fld$benchmark[[i]]
        dt <- data.table(VirusCount=curd$Vir,
                         Mean=me,
                         ReplicateCount=curd$Rep,
                         Variance=curd$Var,
                         curd$sse)
        dff <- rbindlist(list(dff, dt))
      }
      dff
    }))

    df$Model[df$Model == "LMM"] <- "knockdown"
    df$Virus <- paste(df$VirusCount, "viruses")
    df <- dplyr::filter(df, Mean==0.5, Virus=="4 viruses")
    df$Variance[df$Variance == 1] <- "Low variance"
    df$Variance[df$Variance == 2] <- "Medium variance"
    df$Variance[df$Variance == 5] <- "High variance"
    df$Variance <- factor(df$Variance,
                          levels=c("Low variance", "Medium variance" , "High variance"))

    p <-
      ggplot2::ggplot(df) +
      geom_boxplot(aes(x     = Model,
                       y     = MSE,
                       fill  = Model), , width=.5) +
      ggplot2::scale_fill_manual(values =ggthemr:::palettes$fresh$swatch[c(2,4)]) +
      ggplot2::facet_grid(Variance ~ ., scale="free") +
      theme_bw() +
      theme(text = element_text(size = 20)) +
      hrbrthemes::theme_ipsum_rc(base_family="Helvetica") +
      theme(axis.text.x=element_blank(),
            axis.text.y=element_text(size=14),
            axis.title.x=element_blank(),
            axis.title.y=element_text(size=16),
            legend.text=element_text(size=14),
            text = ggplot2::element_text(size = 20),
            strip.text.y = element_text(size = 14)) +
      ylab("Mean residual sum of squares")

    p
}

run <- function()
{
  path <- "./data"
  out.dir <-  "./plots"

  fls <- list.files(path, full.names = T)
  fls.bio <- grep("lmm_predictability_bio", fls, value = T)
  fls.stn <- grep("lmm_predictability_synthetic.*0.5.*", fls, value = T)

  cat("Plotting bio\n")
  p.bio <- bio.pred(fls.bio)
  cat("Plotting syn\n")
  p.syn <- syn.pred(fls.stn)

  leg <- get_legend(p.bio + theme(legend.position="bottom"))
  prow <- plot_grid(p.syn + theme(legend.position="none") +
                      labs(subtitle="Synthetic data benchmark"),
                    p.bio + theme(legend.position="none") +
                      labs(subtitle="Biological data benchmark"),
                    align = 'h',
                    labels = c("(a)", "(b)"),
                    hjust = -1,
                    nrow = 1)
  pl <- plot_grid(prow, leg, ncol=1,  rel_heights= c(10, .5))
  pl

  ggsave(
    filename = paste0(out.dir, "/",
    "predictability_analysis_combined_plot.eps"),
    plot = pl,
    width = 12,
    height = 8)
}

run()
