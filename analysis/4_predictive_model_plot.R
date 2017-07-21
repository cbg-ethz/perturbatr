library(data.table)
library(dtplyr)
library(dplyr)
library(tidyr)
library(lme4)
library(optparse)
library(knockout)
library(ggplot2)
library(uuid)

jaccard <- function(s1, s2) { length(intersect(s1, s2)) / length(union(s1, s2)) }

bio.pred <- function(fls.bio)
{

  stab <- readRDS(fls.bio)
  stab <- stab[-1]

  df <- rbindlist(lapply(stab, function(el) {
    data.table(Virus=length(unlist(stringr::str_split(el$Virus, "_"))),
               el$sse)
  }))

  df <- df %>% dplyr::filter(Sampling=="CV")

  df$Model[df$Model == "LMM"] <- "knockout"
  df$Virus <- paste(df$Virus, "viruses")

  p <-
    ggplot2::ggplot(df) +
    geom_boxplot(aes(x     = Sampling,
                     y     = MSE,
                     fill  = Model)) +
    ggplot2::scale_fill_manual(values = c("#e41a1c", "#377eb8")) +
    ggplot2::facet_grid(Virus ~  .) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 20))


  for (d in dirs)
  {
    ggplot2::ggsave(
      plot=p,
      filename=paste(d, "predictability_bio_data_only_cv.eps", sep = "/"),
      width = 8,
      height = 8
    )
  }

  p

}

syn.pred <- function(fls.stn)
{

    df <- rbindlist(lapply(fls.stn, function (f) {
      fld <- readRDS(f)
      me <- stringr::str_match(f, ".*genemean_(.+)_.*")[2] %>% as.numeric
      dff <- NULL
      # only take lists with benchmark results
      ix <- grep("bt$", names(fld$benchmark))
      for (i in ix)
      {
        print(names(fld$benchmark)[i])
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
        filename=paste(d, "predictability_syn_data_virus_4_replicate_8_genemean_05.eps", sep = "/"),
        plot=p,
        width = 8,
        height = 8
      )
    }

}

###############################################################################


dirs <-
  c("/Users/simondi/PROJECTS/sysvirdrug_project/results/plots",
    "/Users/simondi/PROJECTS/sysvirdrug_project/docs/sysvirdrug_modelling_paper/plots"
  )

fls <- list.files(
  "/Users/simondi/PROJECTS/sysvirdrug_project/results/lmm_predictability_analysis/paper_rds/",
    full.names = T
)

fls.bio <- grep("lmm_predictability_bio",            fls, value = T)
fls.stn <- grep("lmm_predictability__synthetic.*0.5.*", fls, value = T)


bio.pred(fls.bio)
syn.pred(fls.stn)
