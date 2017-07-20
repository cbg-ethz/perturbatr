library(data.table)
library(dtplyr)
library(dplyr)
library(tidyr)
library(lme4)
library(optparse)
library(knockout)
library(ggplot2)
library(uuid)

dirs <-
  c("/Users/simondi/PROJECTS/sysvirdrug_project/results/plots",
    "/Users/simondi/PROJECTS/sysvirdrug_project/docs/sysvirdrug_modelling_paper/plots"
)

fls <-
  list.files(
    "/Users/simondi/PROJECTS/sysvirdrug_project/results/lmm_predictability_analysis/paper_lmm_rds/",
    full.names = T
)
fls.bio <- grep("bio", fls, value = T)
fls.stn <- grep("synthetic", fls, value = T)

bio.pred <- function(fls.bio)
{

  stab <- readRDS(fls.bio)
  stab <- stab[-1]

  df <- rbindlist(lapply(stab, function(el) {
    data.table(Virus=length(unlist(stringr::str_split(el$Virus, "_"))),
               el$sse)
  }))

  df$Model[df$Model == "LMM"] <- "knockout"
  df$Virus <- paste(df$Virus, "viruses")

  p <-
    ggplot2::ggplot(df) +
    geom_boxplot(aes(x     = Sampling,
                     y     = MSE,
                     fill  = Model)) +
    ggplot2::scale_fill_manual(values = c("#e41a1c", "#377eb8")) +
    ggplot2::facet_grid(Virus ~  .) +
    theme_bw() +
    theme(text = element_text(size = 20))
  print(p)

  for (d in dirs)
  {
    ggplot2::ggsave(
      plot=p,
      filename=paste(d, "predictability_bio_data.pdf", sep = "/"),
      width = 8,
      height = 8
    )
  }

}

###############################################################################

jaccard <- function(s1, s2) { length(intersect(s1, s2)) / length(union(s1, s2)) }

syn.pred <- function(fls.stn)
{

  .sse <- function(fls.stn)
  {

    df <- rbindlist(lapply(fls.stn, function (f) {
      fld <- readRDS(f)
      me <- stringr::str_match(f, ".*genemean_(.+)_.*")[2] %>%
        as.numeric
      dff <- NULL
      for (i in seq_along(fld$benchmark))
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

    df$Model[df$Model == "LMM"] <- "knockout"
    df$Virus <- paste(df$VirusCount, "viruses")
    df <- dplyr::filter(df, Mean==0.5, Virus=="4 viruses")
    df$Variance[df$Variance == 1] <- "Low variance"
    df$Variance[df$Variance == 2] <- "Medium variance"
    df$Variance[df$Variance == 5] <- "High variance"
    df$Variance <- factor(df$Variance,
                          levels=c("Low variance", "Medium variance" , "High variance"))
    # df.singlevar <- dplyr::filter(df, Variance==2)
    s <- df
    p <-
      ggplot2::ggplot(df) +
      geom_boxplot(aes(x     = Sampling,
                       y     = MSE,
                       fill  = Model)) +
      ggplot2::scale_fill_manual(values = c("#e41a1c", "#377eb8")) +
      ggplot2::facet_grid(Variance ~  ., scales="free_y") +
      theme_bw() +
      theme(text = element_text(size = 20))
    print(p)

    for (d in dirs)
    {
      ggsave(
        filename=paste(d, "predictability_syn_data_virus_4_replicate_8_genemean_05.pdf", sep = "/"),
        plot=p,
        width = 8,
        height = 8
      )
    }
  }
  .sse(df)

  .jaccard <- function(fls.stn)
  {
    fsvls <- c(10, 25, 50, 75, 100)

    df <- rbindlist(lapply(fls.stn, function (f) {
      fld <- readRDS(f)
      me <- stringr::str_match(f, ".*genemean_(.+)_.*")[2] %>% as.numeric
      dff <- NULL
      for (i in seq_along(fld$benchmark))
      {
        curd <- fld$benchmark[[i]]
        true.dat <- dplyr::select(curd$data, GeneSymbol, GeneEffect) %>% unique %>% .[order(-abs(GeneEffect))]

        btm <- curd$bt.model
        dfff <- NULL
        for (j in seq(length(btm)))
        {
          lme <- btm[[j]]$lmm.fit.effects

          lme.data <- data.table(GeneSymbol = rownames(lme),
                                 GeneEffect = unlist(lme)) %>%
            .[order(-abs(GeneEffect))] %>% as.tbl
          len <- nrow(lme.data)

          m <- data.table(Jaccard=sapply(fsvls, function(i) {
            jaccard(true.dat$GeneSymbol[1:i], lme.data$GeneSymbol[1:i])
          }))
          m$Index <- factor(fsvls, levels=fsvls)
          m$Bootstrap <- j
          dfff <- rbindlist(list(dfff, m))
        }

        dt <- data.table(VirusCount=curd$Vir,
                         Mean=me,
                         ReplicateCount=curd$Rep,
                         Variance=curd$Var,
                         dfff)
        dff <- rbindlist(list(dff, dt))
      }
      dff
    }))
    df$VirusCount <- factor(df$VirusCount, levels=2:4)
    df$Mean <- factor(df$Mean, levels=c(0.5, 1, 2))
    df$Variance <- factor(df$Variance, levels=c(1,2,5))

    l <- df
    df <- l
    df <- dplyr::filter(df, Mean==2)
    ggplot2::ggplot(df, aes(x=Index, y=Jaccard)) +
      ggplot2::geom_boxplot(aes( fill=Variance)) +
      ggplot2::scale_fill_manual(values = c("#e41a1c", "#377eb8", "darkgreen")) +
      ggplot2::facet_grid(VirusCount ~ .) +
      theme_bw() +
      theme(text = element_text(size = 20))
  }


}
