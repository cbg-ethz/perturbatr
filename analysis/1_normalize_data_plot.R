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
library(ggthemr)
library(hrbrthemes)
library(viridis)
library(cowplot)

extrafont::loadfonts()
hrbrthemes::import_roboto_condensed()
ggthemr("fresh", "scientific")

# output directories
dirs <- c("/Users/simondi/PROJECTS/sysvirdrug_project/results/plots",
          "/Users/simondi/PROJECTS/sysvirdrug_project/src/package/analysis/plots",
          "/Users/simondi/PROJECTS/sysvirdrug_project/docs/sysvirdrug_modelling_paper/plots"
)

path <- "/Users/simondi/PROJECTS/sysvirdrug_project/src/package/analysis/"
data.file     <- paste(path, "data/rnai_screen_normalized.rds", sep="/")
x <- readRDS(data.file)

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
  hrbrthemes::theme_ipsum_rc() +
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

for (d in dirs)
{
  ggsave(
    filename = paste(d, "data_overview.eps", sep = "/"),
    plot = pl,
    width = 10,
    height = 10
  )
}
