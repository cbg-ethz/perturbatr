#!/usr/bin/env Rscript

library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(perturbatr)
library(ggplot2)
library(hrbrthemes)
library(cowplot)


ggthemr("fresh", "scientific")


jaccard <- function(s1, s2) { length(intersect(s1, s2)) / length(union(s1, s2)) }


.rank.lmm.benchmark.pair <- function(cur.dat, len.ben, true.ranks, grp, i)
{
  lmm.rel.gen <- true.ranks$GeneSymbol[1:i]
  dff <- NULL
  for (k in 1:(len.ben - 1))
  {
    cur.dat <- tibble::as.tibble(cur.dat)
    e1   <- dplyr::filter(cur.dat, Bootstrap==k) %>%
      dplyr::arrange(desc(abs(Effect)))
    el1  <- dplyr::filter(e1, GeneSymbol %in% lmm.rel.gen) %>%
      dplyr::arrange(GeneSymbol)

    jac1 <- e1[1:i,]$GeneSymbol
    for (l in (k + 1):len.ben)
    {
      e2   <- dplyr::filter(cur.dat, Bootstrap==l) %>%
        dplyr::arrange(desc(abs(Effect)))
      el2  <- dplyr::filter(e2, GeneSymbol %in% lmm.rel.gen) %>%
        dplyr::arrange(GeneSymbol)
      jac2 <- e2[1:i,]$GeneSymbol

      c <- cor(el1$Effect, el2$Effect, method="spearman")
      jac <- jaccard(jac1, jac2)
      dff <- dplyr::bind_rows(list(
        dff,
        as.tibble(rbind(c(grp, GeneCount=i, Score=c,   Method="Spearman"),
                        c(grp, GeneCount=i, Score=jac, Method="Jaccard")))
      ))
    }
  }

  dff
}

rank.lmm.bio <- function(fls.bio, ranking.table, out.dir)
{
  stab            <- readRDS(fls.bio)
  true.ranks      <- stab$full
  bootstrap.data  <- purrr::map_dfr(stab[-1],
    function(e) {
      lmm.genes   <- e[[3]]$lmm.fit$GeneSymbol
      pmm.genes   <- e[[3]]$pmm.fit$GeneSymbol
      lmm.effects <- e[[3]]$lmm.fit$Effect
      pmm.effects <- e[[3]]$pmm.fit$Effect
      tibble(Virus=e[[2]],
             Bootstrap=e[[1]],
             GeneSymbol=c(lmm.genes, pmm.genes),
             Effect=c(lmm.effects, pmm.effects),
             Model=c(rep("perturbatr", length(lmm.genes)),
                     rep("PMM", length(pmm.genes))))
    })

  bootstrap.data <- bootstrap.data %>%
    dplyr::filter(Model=="perturbatr") %>%
    group_by(Virus, Bootstrap) %>%
    dplyr::arrange(-abs(Effect)) %>% ungroup

  virs <- unique(bootstrap.data$Virus)

  true.ranks <- true.ranks$lmm.fit %>%
    dplyr::arrange(desc(abs(Effect)))

  if (!file.exists(ranking.table))
  {
    len.ben <- length(unique(bootstrap.data$Bootstrap))
    df <- NULL
    for (vi in seq_along(virs))
    {
      v <- virs[vi]
      cur.dat <- dplyr::filter(bootstrap.data, Virus == v)
      for (i in c(10, 25, 50, 75, 100))
      {
        dff <- .rank.lmm.benchmark.pair(cur.dat, len.ben, true.ranks, 1, i)
        df  <- bind_rows(
          list(df,
               tibble(Virus=vi,GeneCount=dff$GeneCount,
                 Score=dff$Score,
                 Method=dff$Method)
          )
        )
      }
    }

    saveRDS(df, ranking.table)
  }

  df <- readRDS(ranking.table) %>% as.tibble
  df$Score <- as.numeric(df$Score)
  df$Virus[df$Virus == 1] <- "2 viruses"
  df$Virus[df$Virus == 2] <- "3 viruses"
  df$Virus[df$Virus == 3] <- "4 viruses"

  df$GeneCount <-
    factor(as.integer(as.character(df$GeneCount)), levels=c(10, 25, 50, 75, 100))

  p <-
    ggplot2::ggplot(df) +
    ggplot2::geom_boxplot(ggplot2::aes(GeneCount, Score, fill=Method)) +
    ggplot2::facet_grid(Virus ~ .) +
    ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    hrbrthemes::theme_ipsum_rc(base_family="Helvetica") +
    ggplot2::theme(axis.text.x=element_text(size=14),
          axis.text.y=element_text(size=14),
          axis.title.x=element_text(size=16),
          axis.title.y=element_text(size=16),
          legend.text=element_text(size=14),
          text = ggplot2::element_text(size = 20),
          strip.text.y = element_text(size = 14)) +
    xlab("Number of genes") +
    ylab("Score") +
    guides(fill=guide_legend("Score"))
  p
}


.get.syn.bst.data <- function(fls.stn)
{
  fln.syn.suffixes <- gsub(".*synthetic_data_", "",fls.stn) %>%
    sort() %>%
    unique()
  h <- c()

  # parse the files
  for (f in fls.stn)
  {
    suf <- sub(".*synthetic_data_", "",f)
    h[[suf]] <- f
  }

  h              <- unname(h)
  fl.ls          <- h
  bootstrap.data <- c()

  for (f in fl.ls)
  {
    stab <- readRDS(f)
    stab <- stab[c(-1, -2)]
    stab.data <- c()
    for (i in seq_along(stab))
    {
      el <- stab[[i]]
      lmm.effects <- el$model$lmm.fit$Effect
      pmm.effects <- el$model$pmm.fit$Effect
      lmm.genes   <- el$model$lmm.fit$GeneSymbol
      pmm.genes   <- el$model$pmm.fit$GeneSymbol
      if (el[[2]] >= 7) {
        stab.data <- dplyr::bind_rows(
          list(stab.data,
               tibble(VirusCount=el$Vir,
                      ReplicateCount=el[[2]],
                      Variance=el$Var,
                      Bootstrap=el$Bootstrap,
                      GeneSymbol=c(lmm.genes, pmm.genes),
                      Effect=c(lmm.effects, pmm.effects),
                      Model=c(rep("perturbatr", length(lmm.genes)),
                              rep("PMM", length(pmm.genes))))))
      }
    }

    bootstrap.data <- dplyr::bind_rows(bootstrap.data, stab.data)
  }

  bootstrap.data
}


rank.lmm.syn <- function(fls.stn, ranking.table, out.dir)
{
  bootstrap.data <-
    fls.stn %>%
    .get.syn.bst.data() %>%
    dplyr::filter(Model == "perturbatr") %>%
    dplyr::filter(VirusCount == 4) %>%
    dplyr::group_by(VirusCount, ReplicateCount, Variance) %>%
    { dplyr::mutate(ungroup(.), G = group_indices(.)) }

  if (!file.exists(ranking.table))
  {
    grps <- unique(bootstrap.data$G)
    df <- NULL

    for (grp in grps)
    {
      cur.dat    <- dplyr::filter(bootstrap.data, G == grp)
      true.ranks <- dplyr::filter(cur.dat, Bootstrap == 0) %>%
        dplyr::arrange(desc(abs(Effect)))

      len.ben <- cur.dat$Bootstrap %>%
        unique() %>%
        length()

      for (i in c(10, 25, 50, 75, 100))
      {
        dff <- .rank.lmm.benchmark.pair(cur.dat, len.ben - 1, true.ranks, grp, i)
        df  <- dplyr::bind_rows(
          list(df,
               tibble(
                 VirusCount     = as.integer(cur.dat$VirusCount[1]),
                 ReplicateCount = as.integer(cur.dat$ReplicateCount[1]),
                 Variance       = as.double(cur.dat$Variance[1]),
                 GeneCount      = as.integer(dff$GeneCount),
                 Score          = as.double(dff$Score),
                 Method         = as.character(dff$Method))))
      }
    }

    saveRDS(df, ranking.table)
  }

  df <- readRDS(ranking.table)
  df <- as.tibble(df)

  df$Score <- as.numeric(df$Score)
  df$Variance[df$Variance == 1] <- "Low variance"
  df$Variance[df$Variance == 2] <- "Medium variance"
  df$Variance[df$Variance == 5] <- "High variance"
  df$Variance <- factor(df$Variance,
                        levels=c("Low variance", "Medium variance" , "High variance"))
  df$GeneCount <-
    factor(as.integer(as.character(df$GeneCount)), levels=c(10, 25, 50, 75, 100))

  p <-
    ggplot2::ggplot(df) +
    ggplot2::geom_boxplot(ggplot2::aes(GeneCount, Score, fill=Method)) +
    ggplot2::facet_grid(Variance ~ .) +
    ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size=20)) +
    hrbrthemes::theme_ipsum_rc(base_family="Helvetica") +
    theme(axis.text.x=element_text(size=14),
          axis.text.y=element_text(size=14),
          axis.title.x=element_text(size=16),
          axis.title.y=element_text(size=16),
          legend.text=element_text(size=14),
          text = ggplot2::element_text(size = 20),
          strip.text.y = element_text(size = 14)) +
    xlab("Number of genes") +
    ylab("Score") +
    guides(fill=guide_legend("Score"))

  p
}


run <- function()
{
  data.dir <- "./data"
  out.dir  <- "./plots"
  fls      <- list.files(data.dir, full.names = T)
  fls.bio  <- grep("lmm_stability_bio", fls, value = T)
  fls.stn  <- grep("lmm_stability_synthetic", fls, value = T)

  cat("Plotting bio table\n")
  ranking.table <- paste(data.dir, "lmm_stability_selection_bio_ranking.rds", sep="/")
  p.bio <- rank.lmm.bio(fls.bio, ranking.table, out.dir)

  cat("Plotting syn table\n")
  ranking.table <- paste(data.dir, "lmm_stability_selection_syn_ranking.rds", sep="/")
  p.syn <- rank.lmm.syn(fls.stn, ranking.table, out.dir)

  leg <- get_legend(p.bio + theme(legend.position="bottom"))
  prow <- plot_grid(p.syn + theme(legend.position="none") + labs(subtitle="Synthetic data benchmark"),
                    p.bio + theme(legend.position="none") + labs(subtitle="Biological data benchmark"),
                    align = 'h',
                    labels = c("(a)", "(b)"),
                    hjust = -1,
                    nrow = 1)
  pl <- plot_grid(prow, leg, ncol=1, rel_heights= c(10, .5))

  ggsave(
    filename = paste0(out.dir, "/", "stability_analysis_combined_plot", ".eps"),
    plot = pl,
    width = 14,
    height = 8)
}


run()
