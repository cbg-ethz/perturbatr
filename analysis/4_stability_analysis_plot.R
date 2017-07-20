library(data.table)
library(dtplyr)
library(dplyr)
library(tidyr)
library(lme4)
library(optparse)
library(knockout)
library(ggplot2)
library(uuid)
library(hashmap)

jaccard <- function(s1, s2) { length(intersect(s1, s2)) / length(union(s1, s2)) }

###############################################################################

.rank.lmm.benchmark.pair <- function(cur.dat, len.ben, true.ranks, grp, i)
{
  lmm.rel.gen <- true.ranks$GeneSymbol[1:i]
  dff <- NULL
  for (k in 1:(len.ben-1))
  {
    e1   <- dplyr::filter(cur.dat, Bootstrap==k) %>% .[order(-abs(Effect))]
    el1  <- dplyr::filter(e1, GeneSymbol %in% lmm.rel.gen) %>% .[order(GeneSymbol)]
    jac1 <- e1[1:i]$GeneSymbol
    for (l in (k + 1):len.ben)
    {
      e2   <- dplyr::filter(cur.dat, Bootstrap==l) %>% .[order(-abs(Effect))]
      el2  <- dplyr::filter(e2, GeneSymbol %in% lmm.rel.gen) %>% .[order(GeneSymbol)]
      jac2 <- e2[1:i]$GeneSymbol

      c <- cor(el1$Effect, el2$Effect, method="spearman")
      jac <- jaccard(jac1, jac2)

      dff <- rbindlist(list(
        dff,
        as.data.table(rbind(c(grp , GeneCount=i, Score=c, Method="Spearman"),
                            c(grp , GeneCount=i, Score=jac, Method="Jaccard")))
      ))
    }
  }

  dff
}

.get.syn.bst.data <- function(fls.stn)
{
  # drop some files with same suffices randomly
  # this can be done elegantly using a hash :O
  fln.syn.suffixes <- gsub(".*synthetic_data_", "",fls.stn) %>% sort %>% unique
  h <- c()

  # parse the files
  for (f in fls.stn)
  {
    suf <- sub(".*synthetic_data_", "",f)
    h[[suf]] <- f
  }
  h <- unname(h)

  # in this case we only look at a single file set
  # fl.ls <- grep("viruscnt_6_repcnt_9_var_[1|2|5].rds", h, value=T)
  fl.ls <- h
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
        stab.data <- rbindlist(
          list(stab.data,
               data.table(VirusCount=el$Vir,
                          ReplicateCount=el[[2]],
                          Variance=el$Var,
                          Bootstrap=el$Bootstrap,
                          GeneSymbol=c(lmm.genes, pmm.genes),
                          Effect=c(lmm.effects, pmm.effects),
                          Model=c(rep("knockout", length(lmm.genes)),
                                  rep("PMM", length(pmm.genes)))
               )
          )
        )
      }
    }
    bootstrap.data <- rbindlist(list(bootstrap.data, stab.data))
  }

  bootstrap.data
}

rank.lmm.bio <- function(fls.bio, out.file, ranking.file.prefix)
{

  ranking.table <- paste0(ranking.file.prefix, "_table.rds")
  plt <- paste0(ranking.file.prefix, "_plot.rds")

  stab           <- readRDS(fls.bio)
  true.ranks      <- stab$full
  bootstrap.data <- rbindlist(
    lapply(stab[-1], function(e) {
      lmm.genes   <- e[[3]]$lmm.fit$GeneSymbol
      pmm.genes   <- e[[3]]$pmm.fit$GeneSymbol
      lmm.effects <- e[[3]]$lmm.fit$Effect
      pmm.effects <- e[[3]]$pmm.fit$Effect
      data.table(Virus=e[[2]],
                 Bootstrap=e[[1]],
                 GeneSymbol=c(lmm.genes, pmm.genes),
                 Effect=c(lmm.effects, pmm.effects),
                 Model=c(rep("knockout", length(lmm.genes)),
                         rep("PMM", length(pmm.genes))))
    })
  )

  bootstrap.data <- bootstrap.data %>% dplyr::filter(Model=="knockout")
  bootstrap.data <- bootstrap.data %>%
    group_by(Virus, Bootstrap) %>%
    dplyr::arrange(-abs(Effect)) %>% ungroup
  virs <- unique(bootstrap.data$Virus)

  true.ranks <- true.ranks$lmm.fit %>% .[order(abs(Effect), decreasing=T)]
  if (!file.exists(ranking.table))
  {
    len.ben <- length(unique(bootstrap.data$Bootstrap))
    df <- NULL
    for (vi in seq_along(virs))
    {
      print(paste(vi))
      v <- virs[vi]
      cur.dat <- dplyr::filter(bootstrap.data, Virus==v)
      for (i in c(10, 25, 50, 75, 100))
      {
        print(i)
        dff <- .rank.lmm.benchmark.pair(cur.dat, len.ben, true.ranks, 1, i)
        df  <- rbindlist(
          list(df,
               data.table(
                 Virus=vi,
                 GeneCount=dff$GeneCount,
                 Score=dff$Score,
                 Method=dff$Method)
          )
        )
      }
    }
    saveRDS(df, ranking.table)
  }

  df <- readRDS(ranking.table) %>% as.data.table

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
    ggplot2::scale_fill_manual(values=c("#e41a1c", "#377eb8")) +
    ggplot2::theme_bw() +
    ggplot2::theme(text=element_text(size=20))

  for (d in dirs)
  {
    ggsave(
      filename = paste(d, out.file, sep = "/"),
      plot = p,
      width = 8,
      height = 8
    )
  }

  saveRDS(p, file=plt)

  p
}

rank.lmm.syn <- function(fls.stn, out.file, ranking.file.prefix)
{

  ranking.table <- paste0(ranking.file.prefix, "_table.rds")
  plt <- paste0(ranking.file.prefix, "_plot.rds")

  bootstrap.data <- .get.syn.bst.data(fls.stn)

  bootstrap.data <- bootstrap.data %>% dplyr::filter(Model == "knockout")
  # just look at full data set this time
  bootstrap.data <- dplyr::filter(bootstrap.data, VirusCount==4)
  bootstrap.data <-
    bootstrap.data %>%
    dplyr::group_by(VirusCount, ReplicateCount, Variance) %>%
    dplyr::mutate(G=.GRP)

  if (!file.exists(ranking.table))
  {
    grps <- bootstrap.data$G %>% unique
    df <- NULL
    for (grp in grps)
    {
      cur.dat <- dplyr::filter(bootstrap.data, G == grp)
      true.ranks <- dplyr::filter(cur.dat, Bootstrap == 0) %>%
        dplyr::arrange(-abs(Effect))
      len.ben <- cur.dat$Bootstrap %>% unique %>% length
      for (i in c(10, 25, 50, 75, 100))
      {
        cat(paste(grp, i, "\n"))
        # -1 cause bench goes to 101
        dff <- .rank.lmm.benchmark.pair(cur.dat, len.ben - 1, true.ranks, grp, i)
        df  <- rbindlist(
          list(df,
               data.table(
                 VirusCount=cur.dat$VirusCount[1],
                 ReplicateCount=cur.dat$ReplicateCount[1],
                 Variance=cur.dat$Variance[1],
                 GeneCount=dff$GeneCount,
                 Score=dff$Score,
                 Method=dff$Method)
               )
          )
      }
    }
    saveRDS(df, ranking.table)
  }

  df <- readRDS(ranking.table)
  df <- as.data.table(df)

  df$Score <- as.numeric(df$Score)
  df$Variance[df$Variance == 1] <- "Low variance"
  df$Variance[df$Variance == 2] <- "Medium variance"
  df$Variance[df$Variance == 5] <- "High variance"
  df$Variance <- factor(df$Variance,
                        levels=c("Low variance", "Medium variance" , "High variance"))
  df$GeneCount <-
    factor(as.integer(as.character(df$GeneCount)), levels=c(10, 25, 50, 75, 100))
  # t <- dplyr::filter(df, Virus=="4")

  p <-
    ggplot2::ggplot(df) +
    ggplot2::geom_boxplot(ggplot2::aes(GeneCount, Score, fill=Method)) +
    ggplot2::facet_grid(Variance ~ .) +
    ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::scale_fill_manual(values = c("#e41a1c", "#377eb8")) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size=20))

  for (d in dirs)
  {
    ggsave(
      filename = paste(d, out.file, sep = "/"),
      plot = p,
      width = 8,
      height = 8
    )
  }

  saveRDS(p, file=plt)

  p
}

###############################################################################

# output directories
dirs <- c("/Users/simondi/PROJECTS/sysvirdrug_project/results/plots",
          "/Users/simondi/PROJECTS/sysvirdrug_project/docs/sysvirdrug_modelling_paper/plots"
)

# files to process
fls <-
  list.files(
    "/Users/simondi/PROJECTS/sysvirdrug_project/results/lmm_stability_selection//",
    full.names = T
)

fls.bio <- grep("lmm_stability__bio", fls, value = T)
fls.stn <- grep("lmm_stability__synthetic", fls, value = T)

###############################################################################
# BIOLOGICAL DATA ANALYSIS
#
p <- rank.lmm.bio(fls.bio,
                  "stability_bio_data.eps",
                  "/Users/simondi/PROJECTS/sysvirdrug_project/results/lmm_stability_selection/lmm_stability_selection_bio_ranking")


###############################################################################
# SYNTHETIC DATA ANALYSIS

# rank the stuff
rank.lmm.syn(fls.stn,
             "stability_syn_data.eps",
             "/Users/simondi/PROJECTS/sysvirdrug_project/results/lmm_stability_selection/lmm_stability_selection_syn_ranking")
