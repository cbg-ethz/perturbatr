library(data.table)
library(dtplyr)
library(dplyr)
library(knockout)
library(grid)
library(gridExtra)
library(xtable)
library(igraph)

wd <- "/Users/simondi/PROJECTS/sysvirdrug_project/src/util/knockout_svd_pipeline"

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

network.propagation.single.virus.level <- function(lmm.priorit, graph.file, file.name)
{
  if (!file.exists(graph.file)) graph.file <- "../../data/svd/mappings/fi_flat.tsv"
  gene.path <- dplyr::filter(lmm.priorit$fit$gene.pathogen.effects, FDR <= 0.4)
  gene.called <- do.call("rbind",
                         lapply(unique(gene.path$Virus), function(e) {
                           da <- dplyr::filter(gene.path, Virus==e)
                           diffusion <- knockout::diffuse(da,
                                                          method="mrw",
                                                          path=graph.file,
                                                          r=.8,
                                                          delete.nodes.on.degree=1)$diffusion %>%
                             dplyr::select(GeneSymbol, DiffusionEffect)

                           lmm.hits <- dplyr::select(da, GeneSymbol)
                           genes.not.in.graph <- lmm.hits$GeneSymbol[! da$GeneSymbol %in% diffusion$GeneSymbol]
                           # add htis found by single screen but not in adj matrix
                           diff.hits <- rbindlist(list(diffusion,
                                                       data.table(GeneSymbol=genes.not.in.graph,
                                                                  DiffusionEffect=rep(1, length(genes.not.in.graph)))))

                           data.table(Virus=e, diff.hits)
                         }))

  mrw.summ.hits.single <- dplyr::filter(gene.called, DiffusionEffect == 1)
  mrw.summ.hits.diff   <- dplyr::filter(gene.called, DiffusionEffect != 1) %>%
    dplyr::group_by(Virus) %>%
    dplyr::mutate(Q=quantile(DiffusionEffect, probs=.99)) %>%
    ungroup %>%
    dplyr::filter(DiffusionEffect >= Q) %>%
    dplyr::select(-Q) %>%
    tidyr::spread(Virus, DiffusionEffect) %>%
    as.data.table
  # fill up again
  mrw.summ.hits.diff[is.na(mrw.summ.hits.diff)] <- 0
  mrw.summ.hits.diff <- tidyr::gather(mrw.summ.hits.diff,
                                      Virus, DiffusionEffect,
                                      2:ncol(mrw.summ.hits.diff)) %>%
    as.data.table
  setcolorder(mrw.summ.hits.diff, colnames(mrw.summ.hits.single))

  mrw.fin  <- rbindlist(list(mrw.summ.hits.single, mrw.summ.hits.diff))

  found.count <- mrw.fin %>% dplyr::filter(DiffusionEffect != 0) %>%
    dplyr::group_by(GeneSymbol) %>% dplyr::summarise(n=n())

  found.table <-dplyr::filter(mrw.fin, DiffusionEffect != 0)  %>%
    dplyr::select(-DiffusionEffect) %>%
    tidyr::spread(Virus, Virus) %>% as.data.table

  found.agg <- dplyr::group_by(mrw.summ.hits.diff, GeneSymbol) %>%
    dplyr::summarise(DiffusionEffect=median(DiffusionEffect, na.rm=F))

  diffusion.single.screen.summary <- list(
    Hits=mrw.fin,
    FoundCount= found.count,
    FoundTable=found.table,
    Aggregate=found.agg)

  re.file  <- paste0("hit_selection/single_pathogen_hit_selection/", file.name)
  message(paste("Writing random effect model hits to", re.file))

  saveRDS(diffusion.single.screen.summary, paste0(re.file, ".rds"))

  diffusion.single.screen.summary

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

do.diffusion.on.different.data.sets <- function()
{
  graph.file <-  paste0(wd, "/mappings/fi_flat.tsv")

  primary.lmm <- readRDS(
    paste0(wd, "/hit_selection/all_pathogen_hit_selection/lmm_fit_bootstrap_22_5.rds")
  )

  full.lmm    <- readRDS("hit_selection/all_pathogen_hit_selection/2017_3_24/random_effects_model_hits_full_data.rds")

  primary.diff <-  network.propagation(primary.lmm,
                                       graph.file, "diffusion_primary_data")
  full.diff <-  network.propagation(full.lmm,
                                    graph.file, "diffusion_full_data")

  primary.diff <- readRDS(
    paste0(wd, "/hit_selection/all_pathogen_hit_selection/diffusion_primary_data_remove_1_diff_02_24_5.rds")
  )

  #primary.diff <- readRDS("hit_selection/all_pathogen_hit_selection/diffusion_primary_data_remove_1_diff_02_24_5.rds")
  #full.diff    <- readRDS("hit_selection/all_pathogen_hit_selection/diffusion_full_data.rds")

}

