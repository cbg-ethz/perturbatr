library(data.table)
library(dtplyr)
library(dplyr)
library(knockout)
library(grid)
library(gridExtra)
library(xtable)
library(igraph)

# this is the compilation of the data-sets together
merge.tables <- function()
{

  lmm <- readRDS("~/PROJECTS/sysvirdrug_project/results/hit_selection/all_pathogen_hit_selection/lmm_fit_bootstrap_22_5.rds")
  diffu <- readRDS("~/PROJECTS/sysvirdrug_project/results/hit_selection/all_pathogen_hit_selection/diffusion_primary_data_remove_1_diff_02_24_5.rds")

  lmm.p <- lmm@.gene.effects
  if (all(is.na(lmm.p$Qval))) lmm.p$Qval <- 0.2

  diff.hits <- dplyr::select(diffu@.data, GeneSymbol, DiffusionEffect) %>%
    .[order(-DiffusionEffect)]
  diff.hits <- dplyr::left_join(diff.hits,
                                dplyr::select(lmm.p, GeneSymbol, Qval),
                                by="GeneSymbol")

  gps <- lmm@.gene.pathogen.effects %>%
    dplyr::select(Virus, GeneSymbol, Effect) %>%
    tidyr::spread(Virus, Effect)

  full.table <- diff.hits %>%
    .[order(-abs(DiffusionEffect))] %>%
    dplyr::mutate(RankNetworkAnalysis=seq(.N)) %>%
    dplyr::select(-DiffusionEffect) #%>%
  #.[1:max(which(FDR <= .2))]

  lmm.not.in <-
    lmm.p[!(lmm.p$GeneSymbol %in% full.table$GeneSymbol)] %>%
    dplyr::select(GeneSymbol, Qval) %>% dplyr::mutate(RankNetworkAnalysis = NA_integer_)

  merged.table <- rbindlist(list(lmm.not.in, full.table)) %>%
    dplyr::mutate(o=seq(.N)) %>%
    dplyr::left_join(dplyr::select(lmm@.gene.effects, GeneSymbol, Effect), by="GeneSymbol") %>%
    dplyr::left_join(gps, by="GeneSymbol") %>%
    .[order(o)] %>% dplyr::select(-o)

  setDT(merged.table)[Effect == 0, Effect := NA_real_]




  re.file <- "~/PROJECTS/sysvirdrug_project/results/hit_selection/all_pathogen_hit_selection/diffusion_primary_data_remove_1_diff_01_24_5_merged"
  write.csv(merged.table, paste0(re.file, ".csv"), quote=F, row.names=F)
  merged.table <- read.csv( paste0(re.file, ".csv"), sep=",") %>% as.data.table

  merged.first.25 <-
    merged.table %>%
    dplyr::filter(!is.na(RankNetworkAnalysis)) %>%
    dplyr::rename(Rank=RankNetworkAnalysis) %>%
    dplyr::select(-Qval) %>%
    .[1:25]

  show.xtable(merged.first.25, "caption", "label")

}


sub.graph <- function(diffu)
{
  gene <- "elk1"
  gene <- "dyrk1b"
  gene <- "pkn3"
  gene <- "ssx2ip"

  idx <- which(V(diffu@.graph)$name == gene)
  nei <- as.integer(neighbors(diffu@.graph, gene))
  subgraph <- igraph::induced_subgraph(diffu@.graph, c(idx, nei))
  V(subgraph)$color <- "gray40"
  V(subgraph)[gene]$color <- "darkblue"
  V(subgraph)$size <- 1.25
  V(subgraph)[gene]$size <- 3.25

  plot(subgraph, layout=layout.kamada.kawai, edge.curved = -.05,
       edge.width=2,
       edge.color="grey", vertex.label.family = "Helvetica", vertex.label.font=2, vertex.label.color=V(subgraph)$color,
       vertex.label.cex=V(subgraph)$size, asp=0, rescale = TRUE, vertex.shape="none")
}

network.propagation <- function(lmm, graph.file, file.name)
{
  if (!file.exists(graph.file))
    graph.file <- "../../data/svd/mappings/fi_flat.tsv"

  diffusion <- knockout::diffuse(lmm,
                                 method="mrw",
                                 path=graph.file,
                                 delete.nodes.on.degree=2,
                                 r=0.20)

  re.file  <- paste0(wd, "/hit_selection/all_pathogen_hit_selection/", file.name)
  message(paste("Writing random effect model hits to", re.file))
  saveRDS(diffusion, paste0(re.file, ".rds"))

}

show.diffusion.table <- function(diffusion.hits, label, caption)
{
  best.overall <- diffusion.hits$diffusion %>%
    .[order(DiffusionEffect, decreasing=T)] %>%
    .[1:30] %>%
    as.data.frame
  best.diffusion <- filter(diffusion.hits$diffusion, Effect==0) %>%
    .[order(DiffusionEffect, decreasing=T)] %>%
    .[1:15] %>%
    as.data.frame
  best.overall$Label = "All hits"
  best.diffusion$Label = "Diffusion hits"
  best <- rbind(best.overall, best.diffusion)
  best.list <- split(best, f = best$Label)
  best.list[[1]] <- best.list[[1]][,1:3]
  best.list[[2]] <- best.list[[2]][,1:3]
  attr(best.list, "subheadings") <- paste0(names(best.list))
  bold <- function(x) {paste('{\\textbf{',x,'}}', sep ='')}
  print(xtable(best.overall, caption=caption, label=label),
        digits=5, label=label,
        booktabs=T,
        sanitize.colnames.function=bold,
        sanitize.subheadings.function = bold,
        include.rownames = FALSE,  comment=FALSE)

}


graph.file <-  paste0(wd, "/mappings/fi_flat.tsv")

primary.lmm <- readRDS(
  paste0(wd, "/hit_selection/all_pathogen_hit_selection/lmm_fit_bootstrap_22_5.rds")
)

  full.lmm    <- readRDS(
    paste0(wd, "/hit_selection/all_pathogen_hit_selection/2017_3_24/random_effects_model_hits_full_data.rds")
  )

  primary.diff <-  network.propagation(primary.lmm,
                                       graph.file, "diffusion_primary_data")
  full.diff <-  network.propagation(full.lmm,
                                    graph.file, "diffusion_full_data")

  primary.diff <- readRDS(
    paste0(wd, "/hit_selection/all_pathogen_hit_selection/diffusion_primary_data_remove_1_diff_02_24_5.rds")
  )


