# finish single screen thingins
# do pmm here
# compare pmm to results from single screens (take 'best result') and do this comparison plots again on gene and pathway level
# maybe to the comparison plot of all 6 analyses in one...that would be better. and compare to PMM.
# show most frequent pathways and most frequent genes

library(data.table, quietly=T)
library(dplyr, quietly=T)
library(tidyr, quietly=T)
library(reshape2, quietly=T)
library(svd, quietly=T)

#rm(list = setdiff(ls(), lsf.str()))

load("inst/extdata/rnai.screen.raw.rda")

### hits for chikv using tstatistics and hyperstatistics
# get hits
screen <- (rnai.screen.raw)
chikv.k.dat                 <- filter(screen, Virus== "CHIKV")
sars.k.dat                  <- filter(screen, Virus=="SARS")
hcv.k.dat                   <- filter(screen, Virus=="HCV", Screen=="Kinome")
hcv.g.dat                   <- filter(screen, Virus=="HCV", Screen!="Kinome")
denv.k.dat                  <- filter(screen, Virus=="DENV", Screen=="Kinome")
denv.g.dat                  <- filter(screen, Virus=="DENV", Screen!="Kinome")

class(chikv.k.dat) <- class(sars.k.dat) <- class(hcv.k.dat) <-
  class(hcv.g.dat) <-  class(denv.k.dat) <-  class(denv.g.dat) <- c(class(screen))


# normalize
# TODO normalize SINGLE VIRUS LEVEL
chikv.k.norm <- normalize(chikv.k.dat, normalize=c("robust-z.score"), rm.cytotoxic=T, drop=F)
denv.k.norm  <- normalize(denv.k.dat,  normalize=c("log", "loess", "b.score",  "robust-z.score"), rm.cytotoxic=F, drop=F)

# summarize
screen.summ.loess.bsc.zplate <- summarize(screen.norm.loess.bsc.zplate)
setDT(screen.summ.loess.bsc.zplate)[GeneSymbol == "GAPDH", Control := 0]
setDT(screen.summ.loess.bsc.zplate)[GeneSymbol == "GFP", Control := -0]
setDT(screen.summ.loess.bsc.zplate)[GeneSymbol == "TSG101", Control := 1]
setDT(screen.summ.loess.bsc.zplate)[GeneSymbol == "CD4", Control := -1]
setDT(screen.summ.loess.bsc.zplate)[GeneSymbol == "CD81", Control := 1]
# analyse using PMM
screen.pmms.loess.bsc.zplate <- parallel.mixed.model( screen.summ.loess.bsc.zplate, th=.2)

screen.hits.all.pmm <- as.data.table(screen.pmms.loess.bsc.zplate$res$all.virus.results) %>% .[order(abs(Effect), decreasing=T)]
screen.hits.single.pmm <- do.call("rbind", lapply(screen.pmms.loess.bsc.zplate$res$single.virus.results , function(i) i))
screen.hits.single.pmm$Virus <- unname(unlist(sapply(rownames(screen.hits.single.pmm), function(e) sub(".[[:digit:]]+", "", e))))
single.hit.virus.cpgs <- as.data.table(screen.pmms.loess.bsc.zplate$model$gene.pathogen.effects)

maix <- matrix(0, nrow(screen.hits.all.pmm), length(unique(screen.hits.single.pmm$Virus )))
rownames(maix) <- (screen.hits.all.pmm$GeneSymbol)
colnames(maix) <- (unique(screen.hits.single.pmm$Virus ))
cpg.id <- grep("ccg", colnames(single.hit.virus.cpgs))
row.m <- rownames(maix)
cpg.mat <- filter(single.hit.virus.cpgs, GeneID %in% row.m) %>% dplyr::select(GeneID, grep("ccg", colnames(single.hit.virus.cpgs)))
colnames(cpg.mat) <- c("GeneSymbol", "CHIKV", "DENV", "HCV", "SARS")
class(cpg.mat) <- c(class(cpg.mat),  "svd.cpg.mat")

screen.hits.pmm.max.effect <- screen.hits.all.pmm[order(abs(Effect), decreasing=T)]  %>% .[order(GeneSymbol)]

single.res <- do.call("rbind", lapply(screen.pmms.loess.bsc.zplate$res$single.virus.results, function(i) i))
single.res$Virus <- unname(unlist(sapply(rownames(single.res), function(e) sub(".[[:digit:]]+", "", e))))



chikv.hits.hstat.norm.bsc.zplate              <- prioritize(chikv.hstat.bsc.zplate,   hit.ratio=.6, readout.thresh=2)
sars.hits.hstat.norm.bsc.zplate               <- prioritize(sars.hstat.norm.bsc.zplate, hit.ratio=.6, readout.thresh=2)
hcv.g.hits.hstat.norm.bsc.zplate              <- prioritize(hcv.g.hstat.norm.bsc.zplate, hit.ratio=.6, readout.thresh=2)
hcv.k.hits.hstat.norm.loess.bsc.zplate        <- prioritize(hcv.k.hstat.norm.loess.bsc.zplate, hit.ratio=.6, readout.thresh=1)
denv.k.hits.hstat.norm.loess.bsc.zplate       <- prioritize(denv.k.hstat.norm.loess.bsc.zplate, hit.ratio=.6, readout.thresh=1)
denv.g.hits.hstat.norm.bsc.zplate             <- prioritize(denv.g.hstat.norm.bsc.zplate,     hit.ratio=.6 ,readout.thresh=2)

single.hits.virus.genes  <- c(chikv.hits.hstat.norm.bsc.zplate$Entrez,
                              sars.hits.hstat.norm.bsc.zplate$Entrez,
                              hcv.k.hits.hstat.norm.bsc.zplate$Entrez,
                              hcv.g.hits.hstat.norm.bsc.zplate$Entrez,
                              denv.k.hits.hstat.norm.loess.bsc.zplate$Entrez,
                              denv.g.hits.hstat.norm.bsc.zplate$Entrez)
single.virus.genes <- c(chikv.hstat.bsc.zplate$Entrez,
                        sars.hstat.norm.bsc.zplate$Entrez,
                        hcv.g.hstat.norm.bsc.zplate$Entrez,
                        hcv.k.hstat.norm.loess.bsc.zplate$Entrez,
                        denv.k.hstat.norm.loess.bsc.zplate$Entrez,
                        denv.g.hstat.norm.bsc.zplate$Entrez)

ora.go.single.virus.combined <- ora(single.hits.virus.genes,          single.virus.genes, db="go")
ora.kegg.single.virus.combined <- ora(single.hits.virus.genes,          single.virus.genes, db="kegg")

chikv.ora.single.virus.combined  <- ora(chikv.hits.hstat.norm.bsc.zplate$Entrez,       chikv.hstat.bsc.zplate$Entrez, db="go")$summary
sars.ora.single.virus.combined   <- ora(sars.hits.hstat.norm.bsc.zplate$Entrez,        sars.hstat.norm.bsc.zplate$Entrez, db="go")$summary
hcv.g.ora.single.virus.combined  <- ora(hcv.g.hits.hstat.norm.bsc.zplate$Entrez,       hcv.g.hstat.norm.bsc.zplate$Entrez, db="go")$summary
hcv.k.ora.single.virus.combined  <- ora(hcv.k.hits.hstat.norm.loess.bsc.zplate$Entrez, hcv.k.hstat.norm.loess.bsc.zplate$Entrez, db="go")$summary
denv.g.ora.single.virus.combined <- ora(denv.g.hits.hstat.norm.bsc.zplate$Entrez,      denv.g.hstat.norm.bsc.zplate$Entrez, db="go")$summary
denv.k.ora.single.virus.combined <- ora(denv.k.hits.hstat.norm.loess.bsc.zplate$Entrez,denv.k.hstat.norm.loess.bsc.zplate$Entrez, db="go")$summary

single.virus.genes <- list(CHIKV=chikv.hits.hstat.norm.bsc.zplate$Entrez,
                           DENV_G=denv.g.hits.hstat.norm.bsc.zplate$Entrez,
                           DENV_K=denv.k.hits.hstat.norm.loess.bsc.zplate$Entrez,
                            HCV_G=hcv.g.hits.hstat.norm.bsc.zplate$Entrez,
                            HCV_K=hcv.k.hits.hstat.norm.loess.bsc.zplate$Entrez,
                           SARS=sars.hits.hstat.norm.bsc.zplate$Entrez
                            )

single.vir.ora.terms <- list(CHIKV=chikv.ora.single.virus.combined$Term,
                         SARS=sars.ora.single.virus.combined$Term,
                         HCV_G=hcv.g.ora.single.virus.combined$Term,
                         HCV_K=hcv.k.ora.single.virus.combined$Term,
                         DENV_G=denv.g.ora.single.virus.combined$Term,
                         DENV_K=denv.k.ora.single.virus.combined$Term)
single.vir.go.terms <- list(CHIKV=go.mapping(chikv.hits.hstat.norm.bsc.zplate$Entrez),
                             DENV_G=go.mapping(denv.g.hits.hstat.norm.bsc.zplate$Entrez),
                             DENV_K=go.mapping(denv.k.hits.hstat.norm.loess.bsc.zplate$Entrez),
                             HCV_G=go.mapping(hcv.g.hits.hstat.norm.bsc.zplate$Entrez),
                             HCV_K=go.mapping(hcv.k.hits.hstat.norm.loess.bsc.zplate$Entrez),
                             SARS=go.mapping(sars.hits.hstat.norm.bsc.zplate$Entrez))
single.vir.kegg.terms <- list(CHIKV=kegg.mapping(chikv.hits.hstat.norm.bsc.zplate$Entrez),
                            SARS=kegg.mapping(sars.hits.hstat.norm.bsc.zplate$Entrez),
                            HCV_G=kegg.mapping(hcv.g.hits.hstat.norm.bsc.zplate$Entrez),
                            HCV_K=kegg.mapping(hcv.k.hits.hstat.norm.loess.bsc.zplate$Entrez),
                            DENV_G=kegg.mapping(denv.g.hits.hstat.norm.bsc.zplate$Entrez),
                            DENV_K=kegg.mapping(denv.k.hits.hstat.norm.loess.bsc.zplate$Entrez))

single.vir.gene.overlap <- set.concordance(single.virus.genes)
single.vir.ora.overlap <- set.concordance(single.vir.ora.terms)
single.vir.go.overlap <- set.concordance(single.vir.go.terms)
single.vir.kegg.overlap <- set.concordance(single.vir.kegg.terms)

# all single virus result lists
hit.list <- list(
                  C_K_T_stat_plate       =chikv.hits.tstat.zplate,
                  C_K_T_stat_ctrl        =chikv.hits.tstat.zcont,
                  C_K_T_stat_bsc_plate   =chikv.hits.tstat.bsc.zplate,
                  C_K_T_stat_bsc_ctrl    =chikv.hits.tstat.bsc.zcont,
                  C_K_H_stat_plate       =chikv.hits.hstat.zplate,
                  C_K_H_stat_ctrl        =chikv.hits.hstat.zcont,
                  C_K_H_stat_bsc_plate   =chikv.hits.hstat.bsc.zplate,
                  C_K_H_stat_bsc_ctrl    =chikv.hits.hstat.bsc.zcont,

                  D_G_T_stat_zplate                   =denv.g.hits.tstat.norm.zplate,
                  D_G_T_stat_zctrl                    =denv.g.hits.tstat.norm.zctrl,
                  D_G_T_stat_zctrl_scram              =denv.g.hits.tstat.norm.zctrl.scram,
                  D_G_T_stat_bsc_zplate               =denv.g.hits.tstat.norm.bsc.zplate,
                  D_G_T_stat_bsc_zctrl                =denv.g.hits.tstat.norm.bsc.zctrl,
                  D_G_T_stat_bsc_zplate_scram         =denv.g.hits.tstat.norm.bsc.zctrl.scram,
                  D_G_H_stat_zplate                   =denv.g.hits.hstat.norm.zplate,
                  D_G_H_stat_zctrl                    =denv.g.hits.hstat.norm.zctrl,
                  D_G_H_stat_zctrl_scram              =denv.g.hits.hstat.norm.zctrl.scram,
                  D_G_H_stat_bsc_zplate               =denv.g.hits.hstat.norm.bsc.zplate,
                  D_G_H_stat_bsc_zctrl                =denv.g.hits.hstat.norm.bsc.zctrl,
                  # missed one here
                  D_K_T_stat_zplate                   =denv.k.hits.tstat.norm.zplate,
                  D_K_T_stat_zctrl                    =denv.k.hits.tstat.norm.zctrl,
                  D_K_T_stat_zctrl_scram              =denv.k.hits.tstat.norm.zctrl.scram,
                  D_K_T_stat_bsc_zplate               =denv.k.hits.tstat.norm.bsc.zplate,
                  D_K_T_stat_bsc_zctrl                =denv.k.hits.tstat.norm.bsc.zctrl,
                  D_K_T_stat_bsc_zplate_scram         =denv.k.hits.tstat.norm.bsc.zctrl.scram,
                  D_K_T_stat_loess_bscore_zplate      =denv.k.hits.tstat.norm.loess.bsc.zplate,
                  D_K_T_stat_loess_bscore_zctrl       =denv.k.hits.tstat.norm.loess.bsc.zctrl,
                  D_K_T_stat_loess_bscore_zctrl_scram =denv.k.hits.tstat.norm.loess.bsc.zctrl.scram,
                  D_K_H_stat_zplate                   =denv.k.hits.hstat.norm.zplate,
                  D_K_H_stat_zctrl                    =denv.k.hits.hstat.norm.zctrl,
                  D_K_H_stat_zctrl_scram              =denv.k.hits.hstat.norm.zctrl.scram,
                  D_K_H_stat_bsc_zplate               =denv.k.hits.hstat.norm.bsc.zplate,
                  D_K_H_stat_bsc_zctrl                =denv.k.hits.hstat.norm.bsc.zctrl,
                  D_K_H_stat_bsc_zplate_scram         =denv.k.hits.hstat.norm.bsc.zctrl.scram,
                  D_K_H_stat_loess_bscore_zplate      =denv.k.hits.hstat.norm.loess.bsc.zplate,
                  D_K_H_stat_loess_bscore_zctrl       =denv.k.hits.hstat.norm.loess.bsc.zctrl,
                  D_K_H_stat_loess_bscore_zctrl_scram =denv.k.hits.hstat.norm.loess.bsc.zctrl.scram,

                  H_G_T_stat_zplate                   =hcv.g.hits.tstat.norm.zplate,
                  H_G_T_stat_zctrl                    =hcv.g.hits.tstat.norm.zctrl,
                  H_G_T_stat_zctrl_scram              =hcv.g.hits.tstat.norm.zctrl.scram,
                  H_G_T_stat_bsc_zplate               =hcv.g.hits.tstat.norm.bsc.zplate,
                  H_G_T_stat_bsc_zctrl                =hcv.g.hits.tstat.norm.bsc.zctrl,
                  H_G_T_stat_bsc_zplate_scram         =hcv.g.hits.tstat.norm.bsc.zctrl.scram,
                  H_G_H_stat_zplate                   =hcv.g.hits.hstat.norm.zplate,
                  H_G_H_stat_zctrl                    =hcv.g.hits.hstat.norm.zctrl,
                  H_G_H_stat_zctrl_scram              =hcv.g.hits.hstat.norm.zctrl.scram,
                  H_G_H_stat_bsc_zplate               =hcv.g.hits.hstat.norm.bsc.zplate,
                  H_G_H_stat_bsc_zctrl                =hcv.g.hits.hstat.norm.bsc.zctrl,
                  H_G_H_stat_bsc_zplate_scram         =hcv.g.hits.hstat.norm.bsc.zctrl.scram,

                  H_K_T_stat_zplate                   =hcv.k.hits.tstat.norm.zplate,
                  H_K_T_stat_zctrl                    =hcv.k.hits.tstat.norm.zctrl,
                  H_K_T_stat_zctrl_scram              =hcv.k.hits.tstat.norm.zctrl.scram,
                  H_K_T_stat_bsc_zplate               =hcv.k.hits.tstat.norm.bsc.zplate,
                  H_K_T_stat_bsc_zctrl                =hcv.k.hits.tstat.norm.bsc.zctrl,
                  H_K_T_stat_bsc_zplate_scram         =hcv.k.hits.tstat.norm.bsc.zctrl.scram,
                  H_K_T_stat_loess_bscore_zplate      =hcv.k.hits.tstat.norm.loess.bsc.zplate,
                  H_K_T_stat_loess_bscore_zctrl       =hcv.k.hits.tstat.norm.loess.bsc.zctrl,
                  H_K_T_stat_loess_bscore_zctrl_scram =hcv.k.hits.tstat.norm.loess.bsc.zctrl.scram,
                  H_K_H_stat_zplate                   =hcv.k.hits.hstat.norm.zplate,
                  H_K_H_stat_zctrl                    =hcv.k.hits.hstat.norm.zctrl,
                  H_K_H_stat_zctrl_scram              =hcv.k.hits.hstat.norm.zctrl.scram,
                  H_K_H_stat_bsc_zplate               =hcv.k.hits.hstat.norm.bsc.zplate,
                  H_K_H_stat_bsc_zctrl                =hcv.k.hits.hstat.norm.bsc.zctrl,
                  H_K_H_stat_bsc_zplate_scram         =hcv.k.hits.hstat.norm.bsc.zctrl.scram,
                  H_K_H_stat_loess_bscore_zplate      =hcv.k.hits.hstat.norm.loess.bsc.zplate,
                  H_K_H_stat_loess_bscore_zctrl       =hcv.k.hits.hstat.norm.loess.bsc.zctrl,
                  H_K_H_stat_loess_bscore_zctrl_scram =hcv.k.hits.hstat.norm.loess.bsc.zctrl.scram,

                  S_K_T_stat_zplate     =sars.hits.tstat.norm.zplate,
                  S_K_T_stat_bsc_zplate =sars.hits.tstat.norm.bsc.zplate,
                  S_K_H_stat_zplate     =sars.hits.hstat.norm.zplate,
                  S_K_H_stat_bsc_zplate =sars.hits.hstat.norm.bsc.zplate)


pmm.hit.genes <- data.table(GeneSymbol=screen.pmms.loess.bsc.zplate$res$all.virus.results$GeneSymbol)
gene.ent.map <- dplyr::select(screen.summ.loess.bsc.zplate, GeneSymbol, Entrez) %>% dplyr::filter(!is.na(Entrez)) %>% unique
pmm.hit.genes <- left_join(pmm.hit.genes, gene.ent.map)

pmm.ora.og <- ora(pmm.hit.genes$Entrez, screen.summ.loess.bsc.zplate$Entrez, db="go")
pmm.ora.kegg <- ora(pmm.hit.genes$Entrez, screen.summ.loess.bsc.zplate$Entrez, db="kegg")

#
hit.names <- names(hit.list)
gene.hit.list <- list()
#kegg.hit.list <- list()
for (i in seq(length(hit.list)))
{
  gene.hit.list[[hit.names[i]]] <- hit.list[[i]]$GeneSymbol
#  kegg.na <- kegg.mapping(hit.list[[i]]$Entrez)
 # kegg.hit.list[[hit.names[i]]] <- kegg.na
}

ora.list <- list()
kegg.list <- list()
for (i in seq(length(hit.list)))
{
  ora.go <- ora(unique(as.integer(hit.list[[i]]$Entrez)),
               unique(screen.summ.loess.bsc.zplate$Entrez), db="go")
  ora.kegg <- ora(unique(as.integer(hit.list[[i]]$Entrez)),
               unique(screen.summ.loess.bsc.zplate$Entrez), db="kegg")
  ora.list[[hit.names[i]]] <- ora.go
  kegg.list[[hit.names[i]]] <- ora.kegg
}

saveRDS(ora.list, file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/oras/go_ora.rds")
saveRDS(kegg.list, file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/oras/kegg_ora.rds")

ora.list <- readRDS("/Users/simondi/PHD/data/data/sysvirdrug/target_results/oras/go_ora.rds")
kegg.list <- readRDS("/Users/simondi/PHD/data/data/sysvirdrug/target_results/oras/kegg_ora.rds")

ora.go.gene.list <- list()
ora.go.gene.names <- names(ora.list)
for (i in seq(length(ora.list)))
{
  ora.go.gene.list[[ora.go.gene.names[i]]] <- ora.list[[i]]$summary$Term
}
ora.kegg.gene.list <- list()
ora.kegg.gene.names <- names(ora.list)
for (i in seq(length(kegg.list)))
{
  ora.kegg.gene.list[[ora.kegg.gene.names[i]]] <- kegg.list[[i]]$summary$Term
}

gene.entrez.map <- dplyr::select(screen.summ.loess.bsc.zplate, GeneSymbol, Entrez)
found.genes  <- as.character(screen.pmms.loess.bsc.zplate$res$all.virus.results$GeneSymbol)
found.entrez <- (data.table(GeneSymbol=found.genes) %>%
  dplyr::left_join(gene.entrez.map, by="GeneSymbol") %>%
  unique)$Entrez


all.genes.go.ora   <- ora(unique(as.integer(found.entrez)),
                          unique(screen.summ.loess.bsc.zplate$Entrez), db="go")
all.genes.kegg.ora <- ora(unique(as.integer(found.entrez)),
                          unique(screen.summ.loess.bsc.zplate$Entrez), db="kegg")

jac.gene.hits     <- set.concordance(gene.hit.list)
jac.kegg.hits     <- set.concordance(kegg.hit.list)
jac.kegg.ora.hits <- set.concordance(ora.kegg.gene.list)
jac.go.ora.hits   <- set.concordance(ora.go.gene.list)



df.single.virus.genes <- do.call("rbind",
             lapply(screen.pmms.loess.bsc.zplate$res$single.virus.results,
                    function(i) i))
df.single.virus.genes$Virus <- rownames(df.single.virus.genes)

df.all.virus.genes <- screen.pmms.loess.bsc.zplate$res$all.virus.results

single.res <- do.call("rbind", lapply(screen.pmms.loess.bsc.zplate$res$single.virus.results, function(i) i))
single.res$Virus <- unname(unlist(sapply(rownames(single.res), function(e) sub(".[[:digit:]]+", "", e))))
single.res <-
  group_by(single.res, Virus, Sign=sign(Effect)) %>% dplyr::summarize(cnt=n()) %>% ungroup %>%
  mutate(Count=cnt*Sign)

single.res.box <-
  ggplot(single.res, aes(x=Virus, y=Count, fill=Sign)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_distiller(palette="Spectral") +
  theme_bw() +
  theme(legend.position="none", text = element_text(size=22),
        axis.text.y=element_blank()) +
  geom_hline(yintercept=0) +
  ylab("Count hits") +
  geom_text(aes(x=Virus, y=ifelse(Count>0, Count+1, Count-1), label=abs(Count)), size=5, colour="black")

################ RESULTS PLOTTING

ggsave(single.res.box, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/all_virus_analysis/all_virus_single_hits_bar.pdf")

jac.gene.hit.p1 <- plot(jac.gene.hits, size=20)
ggsave(jac.gene.hit.p1, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/all_virus_analysis/all_virus_analysis_jaccard_gene_hits.pdf")

jac.kegg.hit.p1 <- plot(jac.kegg.hits, size=20)
ggsave(jac.kegg.hit.p1, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/all_virus_analysis/all_virus_analysis_jaccard_kegg_hits.pdf")

jac.kegg.ora.hits.p1 <- plot(jac.kegg.ora.hits, size=6)
ggsave(jac.kegg.ora.hits.p1, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/all_virus_analysis/all_virus_analysis_jaccard_ora_kegg_hits.pdf")

jac.go.ora.hits.p1 <- plot(jac.go.ora.hits, size=6)
ggsave(jac.go.ora.hits.p1, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/all_virus_analysis/all_virus_analysis_jaccard_ora_go_hits.pdf")

# pmm plots
pdf("/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/all_virus_analysis/all_virus_analysis_pmm.pdf")
plot(screen.pmms.loess.bsc.zplate, fdr.thresh=.1, readout.thresh=.05, min.log.scale=F)
dev.off()

write.table(all.genes.go.ora$summary,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/all_virus_analysis/all_virus_analysis_table_combined_gene_list_go_ora.tsv",
            quote=F, row.names=F, sep="\t")


write.table(all.genes.kegg.ora$summary,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/all_virus_analysis/all_virus_analysis_table_combined_gene_list_kegg_ora.tsv",
            quote=F, row.names=F, sep="\t")

write.table(df.single.virus.genes,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/all_virus_analysis/all_virus_analysis_table_single_virus_gene_effects.tsv",
            quote=F, row.names=F, sep="\t")

p <- plot(single.vir.go.overlap, type="o", main="Overlap coefficient on GO-level", size=24)
ggsave(p, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/all_virus_analysis/all_virus_analysis_overlap_coefficient_single_go.pdf",  width=12, height=10)

p <- plot(single.vir.ora.overlap, type="o", main="Overlap coefficient on GO-level after enrichment")
ggsave(p, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/all_virus_analysis/all_virus_analysis_overlap_coefficient_single_go_ora.pdf")

p <- plot(single.vir.gene.overlap, type="o", main="Overlap coefficient on gene-level", size=24 )
ggsave(p, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/all_virus_analysis/all_virus_analysis_overlap_coefficient_single_genes.pdf", width=12, height=10)



data.table(a=rnorm(10),b=rep(0:1, 5))  %>% group_by(b) %>% mutate(g=.GRP)

