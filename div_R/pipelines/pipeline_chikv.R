library(data.table, quietly=T)
library(dplyr, quietly=T)
library(tidyr, quietly=T)
library(reshape2, quietly=T)
load("inst/extdata/rnai.screen.raw.rda")


### hits for chikv using tstatistics and hyperstatistics
# get hits
chikv.dat              <- filter(rnai.screen.raw, Virus=="CHIKV")
# normalize
chikv.norm.zplate       <- normalize(chikv.dat, normalize=c("robust-z.score"), level="plate")
chikv.norm.zcont        <- normalize(chikv.dat, normalize=c("robust-z.score"), level="control")
chikv.norm.bsc.zplate   <- normalize(chikv.dat, normalize=c("b.score", "robust-z.score"), level="plate")
chikv.norm.bsc.zcont    <- normalize(chikv.dat, normalize=c("b.score", "robust-z.score"), level="control")
# summarize
chikv.summ.zplate       <- summarize(chikv.norm.zplate, rm.vial=F)
chikv.summ.zcont        <- summarize(chikv.norm.zcont, rm.vial=F)
chikv.summ.bsc.zplate   <- summarize(chikv.norm.bsc.zplate, rm.vial=F)
chikv.summ.bsc.zcont    <- summarize(chikv.norm.bsc.zcont, rm.vial=F)
# analysis
chikv.tstat.zplate              <- tstatistic(chikv.summ.zplate,      padjust="BH")
chikv.tstat.zcont               <- tstatistic(chikv.summ.zcont,       padjust="BH")
chikv.tstat.bsc.zplate          <- tstatistic(chikv.summ.bsc.zplate,  padjust="BH")
chikv.tstat.bsc.zcont           <- tstatistic(chikv.summ.zcont,       padjust="BH")
chikv.hstat.zplate     <- hyperstatistic(chikv.summ.zplate,           padjust="BH")
chikv.hstat.zcont      <- hyperstatistic(chikv.summ.zcont,            padjust="BH")
chikv.hstat.bsc.zplate <- hyperstatistic(chikv.summ.bsc.zplate,       padjust="BH")
chikv.hstat.bsc.zcont  <- hyperstatistic(chikv.summ.bsc.zcont,        padjust="BH")
# hit calling
chikv.hits.tstat.zplate              <- prioritize(chikv.tstat.zplate,      hit.ratio=.6)
chikv.hits.tstat.zcont               <- prioritize(chikv.tstat.zcont,       hit.ratio=.6)
chikv.hits.tstat.bsc.zplate          <- prioritize(chikv.tstat.bsc.zplate,  hit.ratio=.6)
chikv.hits.tstat.bsc.zcont           <- prioritize(chikv.tstat.bsc.zcont,   hit.ratio=.6)
chikv.hits.hstat.zplate              <- prioritize(chikv.hstat.zplate,      hit.ratio=.6)
chikv.hits.hstat.zcont               <- prioritize(chikv.hstat.zcont,       hit.ratio=.6)
chikv.hits.hstat.bsc.zplate          <- prioritize(chikv.hstat.bsc.zplate,  hit.ratio=.6)
chikv.hits.hstat.bsc.zcont           <- prioritize(chikv.hstat.bsc.zcont,   hit.ratio=.6)

# get number of hits
jac.hits <- set.concordance(list( T_stat_plate       =chikv.hits.tstat.zplate$GeneSymbol,
                                  T_stat_ctrl        =chikv.hits.tstat.zcont$GeneSymbol,
                                  T_stat_bsc_plate   =chikv.hits.tstat.bsc.zplate$GeneSymbol,
                                  T_stat_bsc_ctrl    =chikv.hits.tstat.bsc.zcont$GeneSymbol,
                                  H_stat_plate       =chikv.hits.hstat.zplate$GeneSymbol,
                                  H_stat_ctrl        =chikv.hits.hstat.zcont$GeneSymbol,
                                  H_stat_bsc_plate   =chikv.hits.hstat.bsc.zplate$GeneSymbol,
                                  H_stat_bsc_ctrl    =chikv.hits.hstat.bsc.zcont$GeneSymbol))

jac.kegg.hits <- set.concordance(list( T_stat_plate       =kegg.mapping(chikv.hits.tstat.zplate$Entrez),
                                       T_stat_ctrl        =kegg.mapping(chikv.hits.tstat.zcont$Entrez),
                                       T_stat_bsc_plate   =kegg.mapping(chikv.hits.tstat.bsc.zplate$Entrez),
                                       T_stat_bsc_ctrl    =kegg.mapping(chikv.hits.tstat.bsc.zcont$Entrez),
                                       H_stat_plate       =kegg.mapping(chikv.hits.hstat.zplate$Entrez),
                                       H_stat_ctrl        =kegg.mapping(chikv.hits.hstat.zcont$Entrez),
                                       H_stat_bsc_plate   =kegg.mapping(chikv.hits.hstat.bsc.zplate$Entrez),
                                       H_stat_bsc_ctrl    =kegg.mapping(chikv.hits.hstat.bsc.zcont$Entrez)))

found.genes <- sort(table(c(chikv.hits.tstat.zplate$GeneSymbol,
                            chikv.hits.tstat.zcont$GeneSymbol,
                            chikv.hits.tstat.bsc.zplate$GeneSymbol,
                            chikv.hits.tstat.bsc.zcont$GeneSymbol,
                            chikv.hits.hstat.zplate$GeneSymbol,
                            chikv.hits.hstat.zcont$GeneSymbol,
                            chikv.hits.hstat.bsc.zplate$GeneSymbol,
                            chikv.hits.hstat.bsc.zcont$GeneSymbol)), decreasing=T)

found.entrez <- sort(table(c(chikv.hits.tstat.zplate$Entrez,
                             chikv.hits.tstat.zcont$Entrez,
                             chikv.hits.tstat.bsc.zplate$Entrez,
                             chikv.hits.tstat.bsc.zcont$Entrez,
                             chikv.hits.hstat.zplate$Entrez,
                             chikv.hits.hstat.zcont$Entrez,
                             chikv.hits.hstat.bsc.zplate$Entrez,
                             chikv.hits.hstat.bsc.zcont$Entrez)),
                     decreasing=T)

# ORA
chikv.ora.tstat.zplate               <- ora(chikv.hits.tstat.zplate$Entrez,     chikv.tstat.zplate$Entrez, db="go")
chikv.ora.tstat.zcont                <- ora(chikv.hits.tstat.zcont$Entrez ,     chikv.tstat.zcont$Entrez, db="go")
chikv.ora.tstat.bsc.zplate           <- ora(chikv.hits.tstat.bsc.zplate$Entrez, chikv.tstat.bsc.zplate$Entrez, db="go")
chikv.ora.tstat.bsc.zcont            <- ora(chikv.hits.tstat.bsc.zcont$Entrez,  chikv.tstat.bsc.zcont$Entrez, db="go")
chikv.ora.hstat.zplate               <- ora(chikv.hits.hstat.zplate$Entrez,     chikv.hstat.zplate$Entrez, db="go")
chikv.ora.hstat.zcont                <- ora(chikv.hits.hstat.zcont$Entrez,      chikv.hstat.zcont$Entrez, db="go")
chikv.ora.hstat.bsc.zplate           <- ora(chikv.hits.hstat.bsc.zplate$Entrez, chikv.hstat.bsc.zplate$Entrez, db="go")
chikv.ora.hstat.bsc.zcont            <- ora(chikv.hits.hstat.bsc.zcont$Entrez,  chikv.hstat.bsc.zcont$Entrez, db="go")

all.genes.go.ora  <- ora(unique(as.integer(names(found.entrez))),   chikv.summ.bsc.zcont$Entrez, db="go")
all.genes.kegg.ora <- ora(unique(as.integer(names(found.entrez))),  chikv.summ.bsc.zcont$Entrez, db="kegg")

jac.ora <- set.concordance(list(T_stat_plate    =chikv.ora.tstat.zplate$summary$Term,
                               T_stat_ctrl      =chikv.ora.tstat.zcont$summary$Term,
                               T_stat_bsc_plate =chikv.ora.tstat.bsc.zplate$summary$Term,
                               T_stat_bsc_ctrl  =chikv.ora.tstat.bsc.zcont$summary$Term,
                               H_stat_plate     =chikv.ora.hstat.zplate$summary$Term,
                               H_stat_ctrl      =chikv.ora.hstat.zcont$summary$Term,
                               H_stat_bsc_plate =chikv.ora.hstat.bsc.zplate$summary$Term,
                               H_stat_bsc_ctrl  =chikv.ora.hstat.bsc.zcont$summary$Term))

used.terms <- sort(table(c(chikv.ora.tstat.zplate$summary$Term,
                           chikv.ora.tstat.zcont$summary$Term,
                           chikv.ora.tstat.bsc.zplate$summary$Term,
                           chikv.ora.tstat.bsc.zcont$summary$Term,
                           chikv.ora.hstat.zplate$summary$Term,
                           chikv.ora.hstat.zcont$summary$Term,
                           chikv.ora.hstat.bsc.zplate$summary$Term,
                           chikv.ora.hstat.bsc.zcont$summary$Term)), decreasing=T)


# write all to file
df.used.terms <- data.frame(Number=as.vector(used.terms),   Name=names(used.terms))
df.found.genes <- data.frame(Number=as.vector(found.genes), Name=names(found.genes))

hit.p1 <- plot(jac.hits)
ggsave(hit.p1, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/single_virus_analysis/chikv_kinome_analysis_jaccard_hits.pdf")
ora.p1 <- plot(jac.ora)
ggsave(ora.p1, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/single_virus_analysis/chikv_kinome_analysis_jaccard_ora.pdf")
kegg.p1 <- plot(jac.kegg.hits)
ggsave(kegg.p1, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/single_virus_analysis/chikv_kinome_analysis_jaccard_kegg.pdf")

write.table(df.used.terms,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/chikv_kinome_analysis_table_ora.tsv",
            quote=F, row.names=F, sep="\t")

write.table(df.found.genes,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/chikv_kinome_analysis_table_genes.tsv",
            quote=F, row.names=F, sep="\t")

write.table(all.genes.go.ora$summary,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/chikv_kinome_analysis_table_combined_gene_list_go_ora.tsv",
            quote=F, row.names=F, sep="\t")

write.table(all.genes.kegg.ora$summary,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/chikv_kinome_analysis_table_combined_gene_list_kegg_ora.tsv",
            quote=F, row.names=F, sep="\t")

