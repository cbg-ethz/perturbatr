library(data.table)
library(dtplyr)
library(dplyr)
library(tidyr)
library(lme4)
library(optparse)
library(knockout)
library(ggplot2)
library(uuid)
library(ggthemr)
library(hrbrthemes)
library(viridis)
library(cowplot)

extrafont::loadfonts()
hrbrthemes::import_roboto_condensed()
ggthemr("fresh", "scientific")


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
    geom_boxplot(aes(x     = Model,
                     y     = MSE,
                     fill  = Model), width=.5) +
    ggplot2::facet_grid(Virus ~  .) +
    ggplot2::scale_fill_manual(values =ggthemr:::palettes$fresh$swatch[c(2,4)]) +
    ggplot2::theme_bw() +
    hrbrthemes::theme_ipsum_rc() +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_text(size=14),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=16),
          legend.text=element_text(size=14),
          text = ggplot2::element_text(size = 20),
          strip.text.y = element_text(size = 14)) +
    ylab("Mean residual sum of squares")


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
      geom_boxplot(aes(x     = Model,
                       y     = MSE,
                       fill  = Model), , width=.5) +
      ggplot2::scale_fill_manual(values =ggthemr:::palettes$fresh$swatch[c(2,4)]) +
      ggplot2::facet_grid(Variance ~  ., scale="free") +
      theme_bw() +
      theme(text = element_text(size = 20)) +
      hrbrthemes::theme_ipsum_rc() +
      theme(axis.text.x=element_blank(),
            axis.text.y=element_text(size=14),
            axis.title.x=element_blank(),
            axis.title.y=element_text(size=16),
            legend.text=element_text(size=14),
            text = ggplot2::element_text(size = 20),
            strip.text.y = element_text(size = 14)) +
      ylab("Mean residual sum of squares")

    for (d in dirs)
    {
      ggsave(
        filename=paste(d, "predictability_syn_data_virus_4_replicate_8_genemean_05.eps", sep = "/"),
        plot=p,
        width = 8,
        height = 8
      )
    }

    p

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

fls.bio <- grep("lmm_predictability_bio", fls, value = T)
fls.stn <- grep("lmm_predictability__synthetic.*0.5.*", fls, value = T)


p.bio <- bio.pred(fls.bio)
p.syn <- syn.pred(fls.stn)

leg <- get_legend(p.bio + theme(legend.position="bottom"))
prow <- plot_grid(p.syn + theme(legend.position="none") + labs(subtitle="Synthetic data benchmark"),
                  p.bio + theme(legend.position="none") + labs(subtitle="Biological data benchmark"),
                  align = 'h',
                  labels = c("(a)", "(b)"),
                  hjust = -1,
                  nrow = 1)
pl <- plot_grid(prow, leg, ncol=1,  rel_heights= c(10, .5))
pl

for (d in dirs)
{
  ggsave(
    filename = paste(d, "predictability_analysis_combined_plot.eps", sep = "/"),
    plot = pl,
    width = 12,
    height = 8
  )
}
