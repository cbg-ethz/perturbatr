library(data.table, quietly=T)
library(dplyr, quietly=T)
library(tidyr, quietly=T)
library(reshape2, quietly=T)
load("inst/extdata/rnai.screen.raw.rda")


### hits for chikv using tstatistics and hyperstatistics
# get hits
denv.g.dat              <- filter(rnai.screen.raw, Virus=="DENV", Screen=="DruggableGenome")
# normalize
denv.g.norm.zplate                <- normalize(denv.g.dat, normalize=c("robust-z.score"), level="plate")
denv.g.norm.zctrl                 <- normalize(denv.g.dat, normalize=c("robust-z.score"), level="control")
denv.g.norm.zctrl.scram           <- normalize(denv.g.dat, normalize=c("robust-z.score"), level="control", ctrl="Scrambled")
denv.g.norm.bsc.zplate            <- normalize(denv.g.dat, normalize=c("b.score", "robust-z.score"), level="plate")
denv.g.norm.bsc.zctrl             <- normalize(denv.g.dat, normalize=c("b.score", "robust-z.score"), level="control")
denv.g.norm.bsc.zctrl.scram       <- normalize(denv.g.dat, normalize=c("b.score", "robust-z.score"), level="control", ctrl="Scrambled")
# summarize
denv.g.summ.zplate                <- summarize(denv.g.norm.zplate)
denv.g.summ.zctrl                 <- summarize(denv.g.norm.zctrl)
denv.g.summ.zctrl.scram           <- summarize(denv.g.norm.zctrl.scram)
denv.g.summ.bsc.zplate            <- summarize(denv.g.norm.bsc.zplate)
denv.g.summ.bsc.zctrl             <- summarize(denv.g.norm.bsc.zctrl)
denv.g.summ.bsc.zctrl.scram       <- summarize(denv.g.norm.bsc.zctrl.scram)
# analysis
denv.g.tstat.norm.zplate                   <- tstatistic(denv.g.summ.zplate)
denv.g.tstat.norm.zctrl                    <- tstatistic(denv.g.summ.zctrl,           mu="control")
denv.g.tstat.norm.zctrl.scram              <- tstatistic(denv.g.summ.zctrl.scram,     mu="control", ctrl="Scrambled")
denv.g.tstat.norm.bsc.zplate               <- tstatistic(denv.g.summ.bsc.zplate)
denv.g.tstat.norm.bsc.zctrl                <- tstatistic(denv.g.summ.bsc.zctrl,       mu="control")
denv.g.tstat.norm.bsc.zctrl.scram          <- tstatistic(denv.g.summ.bsc.zctrl.scram, mu="control", ctrl="Scrambled")
denv.g.hstat.norm.zplate                   <- hyperstatistic(denv.g.summ.zplate)
denv.g.hstat.norm.zctrl                    <- hyperstatistic(denv.g.summ.zctrl)
denv.g.hstat.norm.zctrl.scram              <- hyperstatistic(denv.g.summ.zctrl.scram)
denv.g.hstat.norm.bsc.zplate               <- hyperstatistic(denv.g.summ.bsc.zplate)
denv.g.hstat.norm.bsc.zctrl                <- hyperstatistic(denv.g.summ.bsc.zctrl)
denv.g.hstat.norm.bsc.zctrl.scram          <- hyperstatistic(denv.g.summ.bsc.zctrl.scram)
# hit calling
denv.g.hits.tstat.norm.zplate                 <- prioritize(denv.g.tstat.norm.zplate,     hit.ratio=.6)
denv.g.hits.tstat.norm.zctrl                  <- prioritize(denv.g.tstat.norm.zctrl,     hit.ratio=.6)
denv.g.hits.tstat.norm.zctrl.scram            <- prioritize(denv.g.tstat.norm.zctrl.scram,     hit.ratio=.6)
denv.g.hits.tstat.norm.bsc.zplate             <- prioritize(denv.g.tstat.norm.bsc.zplate,     hit.ratio=.6)
denv.g.hits.tstat.norm.bsc.zctrl              <- prioritize(denv.g.tstat.norm.bsc.zctrl,     hit.ratio=.6)
denv.g.hits.tstat.norm.bsc.zctrl.scram        <- prioritize(denv.g.tstat.norm.bsc.zctrl.scram,     hit.ratio=.6)
denv.g.hits.hstat.norm.zplate                 <- prioritize(denv.g.hstat.norm.zplate,     hit.ratio=.6)
denv.g.hits.hstat.norm.zctrl                  <- prioritize(denv.g.hstat.norm.zctrl,     hit.ratio=.6)
denv.g.hits.hstat.norm.zctrl.scram            <- prioritize(denv.g.hstat.norm.zctrl.scram,     hit.ratio=.6)
denv.g.hits.hstat.norm.bsc.zplate             <- prioritize(denv.g.hstat.norm.bsc.zplate,     hit.ratio=.6)
denv.g.hits.hstat.norm.bsc.zctrl              <- prioritize(denv.g.hstat.norm.bsc.zctrl,     hit.ratio=.6)
denv.g.hits.hstat.norm.bsc.zctrl.scram        <- prioritize(denv.g.hstat.norm.bsc.zctrl.scram,     hit.ratio=.6)

# get number of hits
jac.hits <- set.concordance(list( T_stat_zplate                   =denv.g.hits.tstat.norm.zplate$GeneSymbol,
                                  T_stat_zctrl                    =denv.g.hits.tstat.norm.zctrl$GeneSymbol,
                                  T_stat_zctrl_scram              =denv.g.hits.tstat.norm.zctrl.scram$GeneSymbol,
                                  T_stat_bsc_zplate               =denv.g.hits.tstat.norm.bsc.zplate$GeneSymbol,
                                  T_stat_bsc_zctrl                =denv.g.hits.tstat.norm.bsc.zctrl$GeneSymbol,
                                  T_stat_bsc_zplate_scram         =denv.g.hits.tstat.norm.bsc.zctrl.scram$GeneSymbol,
                                  H_stat_zplate                   =denv.g.hits.hstat.norm.zplate$GeneSymbol,
                                  H_stat_zctrl                    =denv.g.hits.hstat.norm.zctrl$GeneSymbol,
                                  H_stat_zctrl_scram              =denv.g.hits.hstat.norm.zctrl.scram$GeneSymbol,
                                  H_stat_bsc_zplate               =denv.g.hits.hstat.norm.bsc.zplate$GeneSymbol,
                                  H_stat_bsc_zctrl                =denv.g.hits.hstat.norm.bsc.zctrl$GeneSymbol,
                                  H_stat_bsc_zplate_scram         =denv.g.hits.hstat.norm.bsc.zctrl.scram$GeneSymbol))
jac.kegg.hits <-  set.concordance(list(T_stat_zplate                   =kegg.mapping(denv.g.hits.tstat.norm.zplate$Entrez),
                                       T_stat_zctrl                    =kegg.mapping(denv.g.hits.tstat.norm.zctrl$Entrez),
                                       T_stat_zctrl_scram              =kegg.mapping(denv.g.hits.tstat.norm.zctrl.scram$Entrez),
                                       T_stat_bsc_zplate               =kegg.mapping(denv.g.hits.tstat.norm.bsc.zplate$Entrez),
                                       T_stat_bsc_zctrl                =kegg.mapping(denv.g.hits.tstat.norm.bsc.zctrl$Entrez),
                                       T_stat_bsc_zplate_scram         =kegg.mapping(denv.g.hits.tstat.norm.bsc.zctrl.scram$Entrez),
                                       H_stat_zplate                   =kegg.mapping(denv.g.hits.hstat.norm.zplate$Entrez),
                                       H_stat_zctrl                    =kegg.mapping(denv.g.hits.hstat.norm.zctrl$Entrez),
                                       H_stat_zctrl_scram              =kegg.mapping(denv.g.hits.hstat.norm.zctrl.scram$Entrez),
                                       H_stat_bsc_zplate               =kegg.mapping(denv.g.hits.hstat.norm.bsc.zplate$Entrez),
                                       H_stat_bsc_zctrl                =kegg.mapping(denv.g.hits.hstat.norm.bsc.zctrl$Entrez),
                                       H_stat_bsc_zplate_scram         =kegg.mapping(denv.g.hits.hstat.norm.bsc.zctrl.scram$Entrez)))
found.genes <- sort(table(c(denv.g.hits.tstat.norm.zplate$GeneSymbol,
                            denv.g.hits.tstat.norm.zctrl$GeneSymbol,
                            denv.g.hits.tstat.norm.zctrl.scram$GeneSymbol,
                            denv.g.hits.tstat.norm.bsc.zplate$GeneSymbol,
                            denv.g.hits.tstat.norm.bsc.zctrl$GeneSymbol,
                            denv.g.hits.tstat.norm.bsc.zctrl.scram$GeneSymbol,
                            denv.g.hits.hstat.norm.zplate$GeneSymbol,
                            denv.g.hits.hstat.norm.zctrl$GeneSymbol,
                            denv.g.hits.hstat.norm.zctrl.scram$GeneSymbol,
                            denv.g.hits.hstat.norm.bsc.zplate$GeneSymbol,
                            denv.g.hits.hstat.norm.bsc.zctrl$GeneSymbol,
                            denv.g.hits.hstat.norm.bsc.zctrl.scram$GeneSymbol)),
                    decreasing=T)
found.entrez <- sort(table(c(denv.g.hits.tstat.norm.zplate$Entrez,
                             denv.g.hits.tstat.norm.zctrl$Entrez,
                             denv.g.hits.tstat.norm.zctrl.scram$Entrez,
                             denv.g.hits.tstat.norm.bsc.zplate$Entrez,
                             denv.g.hits.tstat.norm.bsc.zctrl$Entrez,
                             denv.g.hits.tstat.norm.bsc.zctrl.scram$Entrez,
                             denv.g.hits.hstat.norm.zplate$Entrez,
                             denv.g.hits.hstat.norm.zctrl$Entrez,
                             denv.g.hits.hstat.norm.zctrl.scram$Entrez,
                             denv.g.hits.hstat.norm.bsc.zplate$Entrez,
                             denv.g.hits.hstat.norm.bsc.zctrl$Entrez,
                             denv.g.hits.hstat.norm.bsc.zctrl.scram$Entrez)),
                     decreasing=T)

# ORA
denv.g.ora.tstat.norm.zplate                 <- ora(denv.g.hits.tstat.norm.zplate$Entrez,                denv.g.tstat.norm.zplate$Entrez,                db="go")
denv.g.ora.tstat.norm.zctrl                  <- ora(denv.g.hits.tstat.norm.zctrl$Entrez,                 denv.g.tstat.norm.zctrl$Entrez,                 db="go")
denv.g.ora.tstat.norm.zctrl.scram            <- ora(denv.g.hits.tstat.norm.zctrl.scram$Entrez,           denv.g.tstat.norm.zctrl.scram$Entrez,           db="go")
denv.g.ora.tstat.norm.bsc.zplate             <- ora(denv.g.hits.tstat.norm.bsc.zplate$Entrez,            denv.g.tstat.norm.bsc.zplate$Entrez,            db="go")
denv.g.ora.tstat.norm.bsc.zctrl              <- ora(denv.g.hits.tstat.norm.bsc.zctrl$Entrez,             denv.g.tstat.norm.bsc.zctrl$Entrez,             db="go")
denv.g.ora.tstat.norm.bsc.zctrl.scram        <- ora(denv.g.hits.tstat.norm.bsc.zctrl.scram$Entrez,       denv.g.tstat.norm.bsc.zctrl.scram$Entrez,       db="go")
denv.g.ora.hstat.norm.zplate                 <- ora(denv.g.hits.hstat.norm.zplate$Entrez,                denv.g.hstat.norm.zplate$Entrez,                db="go")
denv.g.ora.hstat.norm.zctrl                  <- ora(denv.g.hits.hstat.norm.zctrl$Entrez,                 denv.g.hstat.norm.zctrl$Entrez,                 db="go")
denv.g.ora.hstat.norm.zctrl.scram            <- ora(denv.g.hits.hstat.norm.zctrl.scram$Entrez,           denv.g.hstat.norm.zctrl.scram$Entrez,           db="go")
denv.g.ora.hstat.norm.bsc.zplate             <- ora(denv.g.hits.hstat.norm.bsc.zplate$Entrez,            denv.g.hstat.norm.bsc.zplate$Entrez,            db="go")
denv.g.ora.hstat.norm.bsc.zctrl              <- ora(denv.g.hits.hstat.norm.bsc.zctrl$Entrez,             denv.g.hstat.norm.bsc.zctrl$Entrez,             db="go")
denv.g.ora.hstat.norm.bsc.zctrl.scram        <- ora(denv.g.hits.hstat.norm.bsc.zctrl.scram$Entrez,       denv.g.hstat.norm.bsc.zctrl.scram$Entrez,       db="go")


all.genes.go.ora   <- ora(unique(as.integer(names(found.entrez))), denv.g.hstat.norm.bsc.zctrl.scram$Entrez, db="go")
all.genes.kegg.ora <- ora(unique(as.integer(names(found.entrez))), denv.g.hstat.norm.bsc.zctrl.scram$Entrez, db="kegg")

jac.ora <- set.concordance(list( T_stat_zplate                    =denv.g.ora.tstat.norm.zplate$summary$Term,
                                 T_stat_zctrl                     =denv.g.ora.tstat.norm.zctrl$summary$Term,
                                 T_stat_zctrl_scram               =denv.g.ora.tstat.norm.zctrl.scram$summary$Term,
                                 T_stat_bsc_zplate                =denv.g.ora.tstat.norm.bsc.zplate$summary$Term,
                                 T_stat_bsc_zctrl                 =denv.g.ora.tstat.norm.bsc.zctrl$summary$Term,
                                 T_stat_bsc_zctrl_scram           =denv.g.ora.tstat.norm.bsc.zctrl.scram$summary$Term,
                                 H_stat_zplate                    =denv.g.ora.hstat.norm.zplate$summary$Term,
                                 H_stat_zctrl                     =denv.g.ora.hstat.norm.zctrl$summary$Term,
                                 H_stat_zctrl_scram               =denv.g.ora.hstat.norm.zctrl.scram$summary$Term,
                                 H_stat_bsc_zplate                =denv.g.ora.hstat.norm.bsc.zplate$summary$Term,
                                 H_stat_bsc_zctrl                 =denv.g.ora.hstat.norm.bsc.zctrl$summary$Term,
                                 H_stat_bsc_zctrl_scram           =denv.g.ora.hstat.norm.bsc.zctrl.scram$summary$Term))

used.terms <- sort(table(c( denv.g.ora.tstat.norm.zplate$summary$Term,
                            denv.g.ora.tstat.norm.zctrl$summary$Term,
                            denv.g.ora.tstat.norm.zctrl.scram$summary$Term,
                            denv.g.ora.tstat.norm.bsc.zplate$summary$Term,
                            denv.g.ora.tstat.norm.bsc.zctrl$summary$Term,
                            denv.g.ora.tstat.norm.bsc.zctrl.scram$summary$Term,
                            denv.g.ora.hstat.norm.zplate$summary$Term,
                            denv.g.ora.hstat.norm.zctrl$summary$Term,
                            denv.g.ora.hstat.norm.zctrl.scram$summary$Term,
                            denv.g.ora.hstat.norm.bsc.zplate$summary$Term,
                            denv.g.ora.hstat.norm.bsc.zctrl$summary$Term,
                            denv.g.ora.hstat.norm.bsc.zctrl.scram$summary$Term)),
                   decreasing=T)

# write all to file
df.used.terms <- data.frame(Number=as.vector(used.terms), Name=names(used.terms))
df.found.genes <- data.frame(Number=as.vector(found.genes), Name=names(found.genes))
hit.p1 <- plot(jac.hits)
ggsave(hit.p1, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/single_virus_analysis/denv_genome_analysis_jaccard_hits.pdf")
ora.p1 <- plot(jac.ora)
ggsave(ora.p1, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/single_virus_analysis/denv_genome_analysis_jaccard_ora.pdf")
kegg.p1 <- plot(jac.kegg.hits)
ggsave(kegg.p1, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/single_virus_analysis/denv_genome_analysis_jaccard_kegg.pdf")

write.table(df.used.terms,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/denv_genome_analysis_table_ora.tsv",
            quote=F, row.names=F, sep="\t")

write.table(df.found.genes,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/denv_genome_analysis_table_genes.tsv",
            quote=F, row.names=F, sep="\t")

write.table(all.genes.go.ora$summary,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/denv_genome_analysis_table_combined_gene_list_go_ora.tsv",
            quote=F, row.names=F, sep="\t")


write.table(all.genes.kegg.ora$summary,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/denv_genome_analysis_table_combined_gene_list_kegg_ora.tsv",
            quote=F, row.names=F, sep="\t")
