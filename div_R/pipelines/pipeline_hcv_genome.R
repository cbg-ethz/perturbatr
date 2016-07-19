library(data.table, quietly=T)
library(dplyr, quietly=T)
library(tidyr, quietly=T)
library(reshape2, quietly=T)

rm(list = setdiff(ls(), lsf.str()))

load("inst/extdata/rnai.screen.raw.rda")

### hits for chikv using tstatistics and hyperstatistics
# get hits
hcv.g.dat              <- filter(rnai.screen.raw, Virus=="HCV", Screen=="DruggableGenome")
# normalize
hcv.g.norm.zplate                <- normalize(hcv.g.dat, normalize=c("robust-z.score"), level="plate")
hcv.g.norm.zctrl                 <- normalize(hcv.g.dat, normalize=c("robust-z.score"), level="control")
hcv.g.norm.zctrl.scram           <- normalize(hcv.g.dat, normalize=c("robust-z.score"), level="control", ctrl="Scrambled")
hcv.g.norm.bsc.zplate            <- normalize(hcv.g.dat, normalize=c("b.score", "robust-z.score"), level="plate")
hcv.g.norm.bsc.zctrl             <- normalize(hcv.g.dat, normalize=c("b.score", "robust-z.score"), level="control")
hcv.g.norm.bsc.zctrl.scram       <- normalize(hcv.g.dat, normalize=c("b.score", "robust-z.score"), level="control", ctrl="Scrambled")
# summarize
hcv.g.summ.zplate                <- summarize(hcv.g.norm.zplate)
hcv.g.summ.zctrl                 <- summarize(hcv.g.norm.zctrl)
hcv.g.summ.zctrl.scram           <- summarize(hcv.g.norm.zctrl.scram)
hcv.g.summ.bsc.zplate            <- summarize(hcv.g.norm.bsc.zplate)
hcv.g.summ.bsc.zctrl             <- summarize(hcv.g.norm.bsc.zctrl)
hcv.g.summ.bsc.zctrl.scram       <- summarize(hcv.g.norm.bsc.zctrl.scram)
# analysis
hcv.g.tstat.norm.zplate                   <- tstatistic(hcv.g.summ.zplate)
hcv.g.tstat.norm.zctrl                    <- tstatistic(hcv.g.summ.zctrl,           mu="control")
hcv.g.tstat.norm.zctrl.scram              <- tstatistic(hcv.g.summ.zctrl.scram,     mu="control", ctrl="Scrambled")
hcv.g.tstat.norm.bsc.zplate               <- tstatistic(hcv.g.summ.bsc.zplate)
hcv.g.tstat.norm.bsc.zctrl                <- tstatistic(hcv.g.summ.bsc.zctrl,       mu="control")
hcv.g.tstat.norm.bsc.zctrl.scram          <- tstatistic(hcv.g.summ.bsc.zctrl.scram, mu="control", ctrl="Scrambled")
hcv.g.hstat.norm.zplate                   <- hyperstatistic(hcv.g.summ.zplate)
hcv.g.hstat.norm.zctrl                    <- hyperstatistic(hcv.g.summ.zctrl)
hcv.g.hstat.norm.zctrl.scram              <- hyperstatistic(hcv.g.summ.zctrl.scram)
hcv.g.hstat.norm.bsc.zplate               <- hyperstatistic(hcv.g.summ.bsc.zplate)
hcv.g.hstat.norm.bsc.zctrl                <- hyperstatistic(hcv.g.summ.bsc.zctrl)
hcv.g.hstat.norm.bsc.zctrl.scram          <- hyperstatistic(hcv.g.summ.bsc.zctrl.scram)
# hit calling
hcv.g.hits.tstat.norm.zplate                 <- prioritize(hcv.g.tstat.norm.zplate,          hit.ratio=.6)
hcv.g.hits.tstat.norm.zctrl                  <- prioritize(hcv.g.tstat.norm.zctrl,           hit.ratio=.6)
hcv.g.hits.tstat.norm.zctrl.scram            <- prioritize(hcv.g.tstat.norm.zctrl.scram,     hit.ratio=.6)
hcv.g.hits.tstat.norm.bsc.zplate             <- prioritize(hcv.g.tstat.norm.bsc.zplate,      hit.ratio=.6)
hcv.g.hits.tstat.norm.bsc.zctrl              <- prioritize(hcv.g.tstat.norm.bsc.zctrl,       hit.ratio=.6)
hcv.g.hits.tstat.norm.bsc.zctrl.scram        <- prioritize(hcv.g.tstat.norm.bsc.zctrl.scram, hit.ratio=.6)
hcv.g.hits.hstat.norm.zplate                 <- prioritize(hcv.g.hstat.norm.zplate,          hit.ratio=.6)
hcv.g.hits.hstat.norm.zctrl                  <- prioritize(hcv.g.hstat.norm.zctrl,           hit.ratio=.6)
hcv.g.hits.hstat.norm.zctrl.scram            <- prioritize(hcv.g.hstat.norm.zctrl.scram,     hit.ratio=.6)
hcv.g.hits.hstat.norm.bsc.zplate             <- prioritize(hcv.g.hstat.norm.bsc.zplate,      hit.ratio=.6)
hcv.g.hits.hstat.norm.bsc.zctrl              <- prioritize(hcv.g.hstat.norm.bsc.zctrl,       hit.ratio=.6)
hcv.g.hits.hstat.norm.bsc.zctrl.scram        <- prioritize(hcv.g.hstat.norm.bsc.zctrl.scram, hit.ratio=.6)

# get number of hits
jac.hits <- set.concordance(list( T_stat_zplate                   =hcv.g.hits.tstat.norm.zplate$GeneSymbol,
                                  T_stat_zctrl                    =hcv.g.hits.tstat.norm.zctrl$GeneSymbol,
                                  T_stat_zctrl_scram              =hcv.g.hits.tstat.norm.zctrl.scram$GeneSymbol,
                                  T_stat_bsc_zplate               =hcv.g.hits.tstat.norm.bsc.zplate$GeneSymbol,
                                  T_stat_bsc_zctrl                =hcv.g.hits.tstat.norm.bsc.zctrl$GeneSymbol,
                                  T_stat_bsc_zplate_scram         =hcv.g.hits.tstat.norm.bsc.zctrl.scram$GeneSymbol,
                                  H_stat_zplate                   =hcv.g.hits.hstat.norm.zplate$GeneSymbol,
                                  H_stat_zctrl                    =hcv.g.hits.hstat.norm.zctrl$GeneSymbol,
                                  H_stat_zctrl_scram              =hcv.g.hits.hstat.norm.zctrl.scram$GeneSymbol,
                                  H_stat_bsc_zplate               =hcv.g.hits.hstat.norm.bsc.zplate$GeneSymbol,
                                  H_stat_bsc_zctrl                =hcv.g.hits.hstat.norm.bsc.zctrl$GeneSymbol,
                                  H_stat_bsc_zplate_scram         =hcv.g.hits.hstat.norm.bsc.zctrl.scram$GeneSymbol))

jac.kegg.hits <-  set.concordance(list(T_stat_zplate                   =kegg.mapping(hcv.g.hits.tstat.norm.zplate$Entrez),
                                       T_stat_zctrl                    =kegg.mapping(hcv.g.hits.tstat.norm.zctrl$Entrez),
                                       T_stat_zctrl_scram              =kegg.mapping(hcv.g.hits.tstat.norm.zctrl.scram$Entrez),
                                       T_stat_bsc_zplate               =kegg.mapping(hcv.g.hits.tstat.norm.bsc.zplate$Entrez),
                                       T_stat_bsc_zctrl                =kegg.mapping(hcv.g.hits.tstat.norm.bsc.zctrl$Entrez),
                                       T_stat_bsc_zplate_scram         =kegg.mapping(hcv.g.hits.tstat.norm.bsc.zctrl.scram$Entrez),
                                       H_stat_zplate                   =kegg.mapping(hcv.g.hits.hstat.norm.zplate$Entrez),
                                       H_stat_zctrl                    =kegg.mapping(hcv.g.hits.hstat.norm.zctrl$Entrez),
                                       H_stat_zctrl_scram              =kegg.mapping(hcv.g.hits.hstat.norm.zctrl.scram$Entrez),
                                       H_stat_bsc_zplate               =kegg.mapping(hcv.g.hits.hstat.norm.bsc.zplate$Entrez),
                                       H_stat_bsc_zctrl                =kegg.mapping(hcv.g.hits.hstat.norm.bsc.zctrl$Entrez),
                                       H_stat_bsc_zplate_scram         =kegg.mapping(hcv.g.hits.hstat.norm.bsc.zctrl.scram$Entrez)))

found.genes <- sort(table(c(hcv.g.hits.tstat.norm.zplate$GeneSymbol,
                            hcv.g.hits.tstat.norm.zctrl$GeneSymbol,
                            hcv.g.hits.tstat.norm.zctrl.scram$GeneSymbol,
                            hcv.g.hits.tstat.norm.bsc.zplate$GeneSymbol,
                            hcv.g.hits.tstat.norm.bsc.zctrl$GeneSymbol,
                            hcv.g.hits.tstat.norm.bsc.zctrl.scram$GeneSymbol,
                            hcv.g.hits.hstat.norm.zplate$GeneSymbol,
                            hcv.g.hits.hstat.norm.zctrl$GeneSymbol,
                            hcv.g.hits.hstat.norm.zctrl.scram$GeneSymbol,
                            hcv.g.hits.hstat.norm.bsc.zplate$GeneSymbol,
                            hcv.g.hits.hstat.norm.bsc.zctrl$GeneSymbol,
                            hcv.g.hits.hstat.norm.bsc.zctrl.scram$GeneSymbol)),
                    decreasing=T)
found.entrez <- sort(table(c(hcv.g.hits.tstat.norm.zplate$Entrez,
                             hcv.g.hits.tstat.norm.zctrl$Entrez,
                             hcv.g.hits.tstat.norm.zctrl.scram$Entrez,
                             hcv.g.hits.tstat.norm.bsc.zplate$Entrez,
                             hcv.g.hits.tstat.norm.bsc.zctrl$Entrez,
                             hcv.g.hits.tstat.norm.bsc.zctrl.scram$Entrez,
                             hcv.g.hits.hstat.norm.zplate$Entrez,
                             hcv.g.hits.hstat.norm.zctrl$Entrez,
                             hcv.g.hits.hstat.norm.zctrl.scram$Entrez,
                             hcv.g.hits.hstat.norm.bsc.zplate$Entrez,
                             hcv.g.hits.hstat.norm.bsc.zctrl$Entrez,
                             hcv.g.hits.hstat.norm.bsc.zctrl.scram$Entrez)),
                     decreasing=T)

# ORA
hcv.g.ora.tstat.norm.zplate                 <- ora(hcv.g.hits.tstat.norm.zplate$Entrez,                hcv.g.tstat.norm.zplate$Entrez,                db="go")
hcv.g.ora.tstat.norm.zctrl                  <- ora(hcv.g.hits.tstat.norm.zctrl$Entrez,                 hcv.g.tstat.norm.zctrl$Entrez,                 db="go")
hcv.g.ora.tstat.norm.zctrl.scram            <- ora(hcv.g.hits.tstat.norm.zctrl.scram$Entrez,           hcv.g.tstat.norm.zctrl.scram$Entrez,           db="go")
hcv.g.ora.tstat.norm.bsc.zplate             <- ora(hcv.g.hits.tstat.norm.bsc.zplate$Entrez,            hcv.g.tstat.norm.bsc.zplate$Entrez,            db="go")
hcv.g.ora.tstat.norm.bsc.zctrl              <- ora(hcv.g.hits.tstat.norm.bsc.zctrl$Entrez,             hcv.g.tstat.norm.bsc.zctrl$Entrez,             db="go")
hcv.g.ora.tstat.norm.bsc.zctrl.scram        <- ora(hcv.g.hits.tstat.norm.bsc.zctrl.scram$Entrez,       hcv.g.tstat.norm.bsc.zctrl.scram$Entrez,       db="go")
hcv.g.ora.hstat.norm.zplate                 <- ora(hcv.g.hits.hstat.norm.zplate$Entrez,                hcv.g.hstat.norm.zplate$Entrez,                db="go")
hcv.g.ora.hstat.norm.zctrl                  <- ora(hcv.g.hits.hstat.norm.zctrl$Entrez,                 hcv.g.hstat.norm.zctrl$Entrez,                 db="go")
hcv.g.ora.hstat.norm.zctrl.scram            <- ora(hcv.g.hits.hstat.norm.zctrl.scram$Entrez,           hcv.g.hstat.norm.zctrl.scram$Entrez,           db="go")
hcv.g.ora.hstat.norm.bsc.zplate             <- ora(hcv.g.hits.hstat.norm.bsc.zplate$Entrez,            hcv.g.hstat.norm.bsc.zplate$Entrez,            db="go")
hcv.g.ora.hstat.norm.bsc.zctrl              <- ora(hcv.g.hits.hstat.norm.bsc.zctrl$Entrez,             hcv.g.hstat.norm.bsc.zctrl$Entrez,             db="go")
hcv.g.ora.hstat.norm.bsc.zctrl.scram        <- ora(hcv.g.hits.hstat.norm.bsc.zctrl.scram$Entrez,       hcv.g.hstat.norm.bsc.zctrl.scram$Entrez,       db="go")


all.genes.go.ora   <- ora(unique(as.integer(names(found.entrez))), hcv.g.hstat.norm.bsc.zctrl.scram$Entrez, db="go")
all.genes.kegg.ora <- ora(unique(as.integer(names(found.entrez))), hcv.g.hstat.norm.bsc.zctrl.scram$Entrez, db="kegg")

jac.ora <- set.concordance(list( T_stat_zplate                    =hcv.g.ora.tstat.norm.zplate$summary$Term,
                                 T_stat_zctrl                     =hcv.g.ora.tstat.norm.zctrl$summary$Term,
                                 T_stat_zctrl_scram               =hcv.g.ora.tstat.norm.zctrl.scram$summary$Term,
                                 T_stat_bsc_zplate                =hcv.g.ora.tstat.norm.bsc.zplate$summary$Term,
                                 T_stat_bsc_zctrl                 =hcv.g.ora.tstat.norm.bsc.zctrl$summary$Term,
                                 T_stat_bsc_zctrl_scram           =hcv.g.ora.tstat.norm.bsc.zctrl.scram$summary$Term,
                                 H_stat_zplate                    =hcv.g.ora.hstat.norm.zplate$summary$Term,
                                 H_stat_zctrl                     =hcv.g.ora.hstat.norm.zctrl$summary$Term,
                                 H_stat_zctrl_scram               =hcv.g.ora.hstat.norm.zctrl.scram$summary$Term,
                                 H_stat_bsc_zplate                =hcv.g.ora.hstat.norm.bsc.zplate$summary$Term,
                                 H_stat_bsc_zctrl                 =hcv.g.ora.hstat.norm.bsc.zctrl$summary$Term,
                                 H_stat_bsc_zctrl_scram           =hcv.g.ora.hstat.norm.bsc.zctrl.scram$summary$Term))

used.terms <- sort(table(c( hcv.g.ora.tstat.norm.zplate$summary$Term,
                            hcv.g.ora.tstat.norm.zctrl$summary$Term,
                            hcv.g.ora.tstat.norm.zctrl.scram$summary$Term,
                            hcv.g.ora.tstat.norm.bsc.zplate$summary$Term,
                            hcv.g.ora.tstat.norm.bsc.zctrl$summary$Term,
                            hcv.g.ora.tstat.norm.bsc.zctrl.scram$summary$Term,
                            hcv.g.ora.hstat.norm.zplate$summary$Term,
                            hcv.g.ora.hstat.norm.zctrl$summary$Term,
                            hcv.g.ora.hstat.norm.zctrl.scram$summary$Term,
                            hcv.g.ora.hstat.norm.bsc.zplate$summary$Term,
                            hcv.g.ora.hstat.norm.bsc.zctrl$summary$Term,
                            hcv.g.ora.hstat.norm.bsc.zctrl.scram$summary$Term)),
                   decreasing=T)

# write all to file
df.used.terms <- data.frame(Number=as.vector(used.terms), Name=names(used.terms))
df.found.genes <- data.frame(Number=as.vector(found.genes), Name=names(found.genes))
hit.p1 <- plot(jac.hits)
ggsave(hit.p1, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/single_virus_analysis/hcv_genome_analysis_jaccard_hits.pdf")
ora.p1 <- plot(jac.ora)
ggsave(ora.p1, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/single_virus_analysis/hcv_genome_analysis_jaccard_ora.pdf")
kegg.p1 <- plot(jac.kegg.hits)
ggsave(kegg.p1, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/single_virus_analysis/hcv_genome_analysis_jaccard_kegg.pdf")

write.table(df.used.terms,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/hcv_genome_analysis_table_ora.tsv",
            quote=F, row.names=F, sep="\t")

write.table(df.found.genes,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/hcv_genome_analysis_table_genes.tsv",
            quote=F, row.names=F, sep="\t")

write.table(all.genes.go.ora$summary,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/hcv_genome_analysis_table_combined_gene_list_go_ora.tsv",
            quote=F, row.names=F, sep="\t")

write.table(all.genes.kegg.ora$summary,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/hcv_genome_analysis_table_combined_gene_list_kegg_ora.tsv",
            quote=F, row.names=F, sep="\t")
