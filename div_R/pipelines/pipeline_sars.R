library(data.table, quietly=T)
library(dplyr, quietly=T)
library(tidyr, quietly=T)
library(reshape2, quietly=T)
load("inst/extdata/rnai.screen.raw.rda")


### hits for chikv using tstatistics and hyperstatistics
# get hits
sars.dat              <- filter(rnai.screen.raw, Virus=="SARS")
# normalize
sars.norm.zplate       <- normalize(sars.dat, normalize=c("robust-z.score"), level="plate")
sars.norm.bsc.zplate   <- normalize(sars.dat, normalize=c("b.score", "robust-z.score"), level="plate")
# summarize
sars.summ.zplate       <- summarize(sars.norm.zplate,     rm.vial=F)
sars.summ.bsc.zplate   <- summarize(sars.norm.bsc.zplate, rm.vial=F)
# analysis
sars.tstat.norm.zplate         <- tstatistic(    sars.summ.zplate,     padjust="BH")
sars.tstat.norm.bsc.zplate     <- tstatistic(    sars.summ.bsc.zplate, padjust="BH")
sars.hstat.norm.zplate         <- hyperstatistic(sars.summ.zplate,     padjust="BH")
sars.hstat.norm.bsc.zplate     <- hyperstatistic(sars.summ.bsc.zplate, padjust="BH")
# hit calling
sars.hits.tstat.norm.zplate     <- prioritize(sars.tstat.norm.zplate,     hit.ratio=.6)
sars.hits.tstat.norm.bsc.zplate <- prioritize(sars.tstat.norm.bsc.zplate, hit.ratio=.6)
sars.hits.hstat.norm.zplate     <- prioritize(sars.hstat.norm.zplate,     hit.ratio=.6)
sars.hits.hstat.norm.bsc.zplate <- prioritize(sars.hstat.norm.bsc.zplate, hit.ratio=.6)

# get number of hits
jac.hits <- set.concordance(list( T_stat_zplate     =sars.hits.tstat.norm.zplate$GeneSymbol,
                                  T_stat_bsc_zplate =sars.hits.tstat.norm.bsc.zplate$GeneSymbol,
                                  H_stat_zplate     =sars.hits.hstat.norm.zplate$GeneSymbol,
                                  H_stat_bsc_zplate =sars.hits.hstat.norm.bsc.zplate$GeneSymbol))

jac.kegg.hits <-  set.concordance(list(T_stat_zplate                   =kegg.mapping(sars.hits.tstat.norm.zplate$Entrez),
                                       T_stat_bsc_zplate               =kegg.mapping(sars.hits.tstat.norm.bsc.zplate$Entrez),
                                       H_stat_zplate                   =kegg.mapping(sars.hits.hstat.norm.zplate$Entrez),
                                       H_stat_bsc_zplate               =kegg.mapping(sars.hits.hstat.norm.bsc.zplate$Entrez)))

found.genes <- sort(table(c(sars.hits.tstat.norm.zplate$GeneSymbol,
                            sars.hits.tstat.norm.bsc.zplate$GeneSymbol,
                            sars.hits.hstat.norm.zplate$GeneSymbol,
                            sars.hits.hstat.norm.bsc.zplate$GeneSymbol)), decreasing=T)

found.entrez <- sort(table(c(sars.hits.tstat.norm.zplate$Entrez,
                             sars.hits.tstat.norm.bsc.zplate$Entrez,
                             sars.hits.hstat.norm.zplate$Entrez,
                             sars.hits.hstat.norm.bsc.zplate$Entrez)), decreasing=T)

# ora
sars.ora.tstat.zplate      <- ora(sars.hits.tstat.norm.zplate$Entrez,     sars.tstat.norm.zplate$Entrez,     db="go")
sars.ora.tstat.bsc.zplate  <- ora(sars.hits.tstat.norm.bsc.zplate$Entrez, sars.tstat.norm.bsc.zplate$Entrez, db="go")
sars.ora.hstat.zplate      <- ora(sars.hits.hstat.norm.zplate$Entrez,     sars.hstat.norm.zplate$Entrez,     db="go")
sars.ora.hstat.bsc.zplate  <- ora(sars.hits.hstat.norm.bsc.zplate$Entrez, sars.hstat.norm.bsc.zplate$Entrez, db="go")

all.genes.go.ora   <- ora(unique(as.integer(names(found.entrez))), sars.hstat.norm.zplate$Entrez, db="go")
all.genes.kegg.ora <- ora(unique(as.integer(names(found.entrez))), sars.hstat.norm.zplate$Entrez, db="kegg")


jac.ora <- set.concordance(list( T_stat_zplate     =sars.ora.tstat.zplate$summary$Term,
                                 T_stat_bsc_zplate =sars.ora.tstat.bsc.zplate$summary$Term,
                                 H_stat_zplate     =sars.ora.hstat.zplate$summary$Term,
                                 H_stat_bsc_zplate =sars.ora.hstat.bsc.zplate$summary$Term))

used.terms <- sort(table(c(sars.ora.tstat.zplate$summary$Term,
                           sars.ora.tstat.bsc.zplate$summary$Term,
                           sars.ora.hstat.zplate$summary$Term,
                           sars.ora.hstat.bsc.zplate$summary$Term)),
                   decreasing=T)


# write all to file
df.used.terms <- data.frame(Number=as.vector(used.terms), Name=names(used.terms))
df.found.genes <- data.frame(Number=as.vector(found.genes), Name=names(found.genes))

hit.p1 <- plot(jac.hits)
ggsave(hit.p1, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/single_virus_analysis/sars_kinome_analysis_jaccard_hits.pdf")

ora.p1 <- plot(jac.ora)
ggsave(ora.p1, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/single_virus_analysis/sars_kinome_analysis_jaccard_ora.pdf")

kegg.p1 <- plot(jac.kegg.hits)
ggsave(kegg.p1, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/single_virus_analysis/sars_kinome_analysis_jaccard_kegg.pdf")

write.table(df.used.terms,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/sars_kinome_analysis_table_ora.tsv",
            quote=F, row.names=F, sep="\t")

write.table(df.found.genes,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/sars_kinome_analysis_table_genes.tsv",
            quote=F, row.names=F, sep="\t")

write.table(all.genes.go.ora$summary,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/sars_kinome_analysis_table_combined_gene_list_go_ora.tsv",
            quote=F, row.names=F, sep="\t")


write.table(all.genes.kegg.ora$summary,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/sars_kinome_analysis_table_combined_gene_list_kegg_ora.tsv",
            quote=F, row.names=F, sep="\t")
