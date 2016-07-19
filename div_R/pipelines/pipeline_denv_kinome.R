library(data.table, quietly=T)
library(dplyr, quietly=T)
library(tidyr, quietly=T)
library(reshape2, quietly=T)
load("inst/extdata/rnai.screen.raw.rda")

### hits for chikv using tstatistics and hyperstatistics
# get hits
denv.k.dat              <- filter(rnai.screen.raw, Virus=="DENV", Screen=="Kinome")
# normalize
denv.k.norm.zplate                <- normalize(denv.k.dat, normalize=c("robust-z.score"), level="plate")
denv.k.norm.zctrl                 <- normalize(denv.k.dat, normalize=c("robust-z.score"), level="control")
denv.k.norm.zctrl.scram           <- normalize(denv.k.dat, normalize=c("robust-z.score"), level="control", ctrl="Scrambled")
denv.k.norm.bsc.zplate            <- normalize(denv.k.dat, normalize=c("b.score", "robust-z.score"), level="plate")
denv.k.norm.bsc.zctrl             <- normalize(denv.k.dat, normalize=c("b.score", "robust-z.score"), level="control")
denv.k.norm.bsc.zctrl.scram       <- normalize(denv.k.dat, normalize=c("b.score", "robust-z.score"), level="control", ctrl="Scrambled")
denv.k.norm.loess.bsc.zplate      <- normalize(denv.k.dat, normalize=c("loess", "b.score", "robust-z.score"), level="plate")
denv.k.norm.loess.bsc.zctrl       <- normalize(denv.k.dat, normalize=c("loess", "b.score", "robust-z.score"), level="control")
denv.k.norm.loess.bsc.zctrl.scram <- normalize(denv.k.dat, normalize=c("loess", "b.score", "robust-z.score"), level="control", ctrl="Scrambled")
# summarize
denv.k.summ.zplate                <- summarize(denv.k.norm.zplate)
denv.k.summ.zctrl                 <- summarize(denv.k.norm.zctrl)
denv.k.summ.zctrl.scram           <- summarize(denv.k.norm.zctrl.scram)
denv.k.summ.bsc.zplate            <- summarize(denv.k.norm.bsc.zplate)
denv.k.summ.bsc.zctrl             <- summarize(denv.k.norm.bsc.zctrl)
denv.k.summ.bsc.zctrl.scram       <- summarize(denv.k.norm.bsc.zctrl.scram)
denv.k.summ.loess.bsc.zplate      <- summarize(denv.k.norm.loess.bsc.zplate)
denv.k.summ.loess.bsc.zctrl       <- summarize(denv.k.norm.loess.bsc.zctrl)
denv.k.summ.loess.bsc.zctrl.scram <- summarize(denv.k.norm.loess.bsc.zctrl.scram)
# analysis
denv.k.tstat.norm.zplate                   <- tstatistic(denv.k.summ.zplate,                padjust="BH")
denv.k.tstat.norm.zctrl                    <- tstatistic(denv.k.summ.zctrl,                 padjust="BH", mu="control")
denv.k.tstat.norm.zctrl.scram              <- tstatistic(denv.k.summ.zctrl.scram,           padjust="BH", mu="control", ctrl="Scrambled")
denv.k.tstat.norm.bsc.zplate               <- tstatistic(denv.k.summ.bsc.zplate,            padjust="BH")
denv.k.tstat.norm.bsc.zctrl                <- tstatistic(denv.k.summ.bsc.zctrl,             padjust="BH", mu="control")
denv.k.tstat.norm.bsc.zctrl.scram          <- tstatistic(denv.k.summ.bsc.zctrl.scram,       padjust="BH", mu="control", ctrl="Scrambled")
denv.k.tstat.norm.loess.bsc.zplate         <- tstatistic(denv.k.summ.loess.bsc.zplate,      padjust="BH")
denv.k.tstat.norm.loess.bsc.zctrl          <- tstatistic(denv.k.summ.loess.bsc.zctrl,       padjust="BH", mu="control")
denv.k.tstat.norm.loess.bsc.zctrl.scram    <- tstatistic(denv.k.summ.loess.bsc.zctrl.scram, padjust="BH", mu="control", ctrl="Scrambled")
denv.k.hstat.norm.zplate                   <- hyperstatistic(denv.k.summ.zplate)
denv.k.hstat.norm.zctrl                    <- hyperstatistic(denv.k.summ.zctrl)
denv.k.hstat.norm.zctrl.scram              <- hyperstatistic(denv.k.summ.zctrl.scram)
denv.k.hstat.norm.bsc.zplate               <- hyperstatistic(denv.k.summ.bsc.zplate)
denv.k.hstat.norm.bsc.zctrl                <- hyperstatistic(denv.k.summ.bsc.zctrl)
denv.k.hstat.norm.bsc.zctrl.scram          <- hyperstatistic(denv.k.summ.bsc.zctrl.scram)
denv.k.hstat.norm.loess.bsc.zplate         <- hyperstatistic(denv.k.summ.loess.bsc.zplate)
denv.k.hstat.norm.loess.bsc.zctrl          <- hyperstatistic(denv.k.summ.loess.bsc.zctrl)
denv.k.hstat.norm.loess.bsc.zctrl.scram    <- hyperstatistic(denv.k.summ.loess.bsc.zctrl.scram)
# hit calling
denv.k.hits.tstat.norm.zplate                 <- prioritize(denv.k.tstat.norm.zplate,                hit.ratio=.6)
denv.k.hits.tstat.norm.zctrl                  <- prioritize(denv.k.tstat.norm.zctrl,                 hit.ratio=.6)
denv.k.hits.tstat.norm.zctrl.scram            <- prioritize(denv.k.tstat.norm.zctrl.scram,           hit.ratio=.6)
denv.k.hits.tstat.norm.bsc.zplate             <- prioritize(denv.k.tstat.norm.bsc.zplate,            hit.ratio=.6)
denv.k.hits.tstat.norm.bsc.zctrl              <- prioritize(denv.k.tstat.norm.bsc.zctrl,             hit.ratio=.6)
denv.k.hits.tstat.norm.bsc.zctrl.scram        <- prioritize(denv.k.tstat.norm.bsc.zctrl.scram,       hit.ratio=.6)
denv.k.hits.tstat.norm.loess.bsc.zplate       <- prioritize(denv.k.tstat.norm.loess.bsc.zplate,      hit.ratio=.6)
denv.k.hits.tstat.norm.loess.bsc.zctrl        <- prioritize(denv.k.tstat.norm.loess.bsc.zctrl,       hit.ratio=.6)
denv.k.hits.tstat.norm.loess.bsc.zctrl.scram  <- prioritize(denv.k.tstat.norm.loess.bsc.zctrl.scram, hit.ratio=.6)
denv.k.hits.hstat.norm.zplate                 <- prioritize(denv.k.hstat.norm.zplate,                hit.ratio=.6)
denv.k.hits.hstat.norm.zctrl                  <- prioritize(denv.k.hstat.norm.zctrl,                 hit.ratio=.6)
denv.k.hits.hstat.norm.zctrl.scram            <- prioritize(denv.k.hstat.norm.zctrl.scram,           hit.ratio=.6)
denv.k.hits.hstat.norm.bsc.zplate             <- prioritize(denv.k.hstat.norm.bsc.zplate,            hit.ratio=.6)
denv.k.hits.hstat.norm.bsc.zctrl              <- prioritize(denv.k.hstat.norm.bsc.zctrl,             hit.ratio=.6)
denv.k.hits.hstat.norm.bsc.zctrl.scram        <- prioritize(denv.k.hstat.norm.bsc.zctrl.scram,       hit.ratio=.6)
denv.k.hits.hstat.norm.loess.bsc.zplate       <- prioritize(denv.k.hstat.norm.loess.bsc.zplate,      hit.ratio=.6)
denv.k.hits.hstat.norm.loess.bsc.zctrl        <- prioritize(denv.k.hstat.norm.loess.bsc.zctrl,       hit.ratio=.6)
denv.k.hits.hstat.norm.loess.bsc.zctrl.scram  <- prioritize(denv.k.hstat.norm.loess.bsc.zctrl.scram, hit.ratio=.6)
# get number of hits
jac.hits <- set.concordance(list( T_stat_zplate                   =denv.k.hits.tstat.norm.zplate$GeneSymbol,
                                  T_stat_zctrl                    =denv.k.hits.tstat.norm.zctrl$GeneSymbol,
                                  T_stat_zctrl_scram              =denv.k.hits.tstat.norm.zctrl.scram$GeneSymbol,
                                  T_stat_bsc_zplate               =denv.k.hits.tstat.norm.bsc.zplate$GeneSymbol,
                                  T_stat_bsc_zctrl                =denv.k.hits.tstat.norm.bsc.zctrl$GeneSymbol,
                                  T_stat_bsc_zplate_scram         =denv.k.hits.tstat.norm.bsc.zctrl.scram$GeneSymbol,
                                  T_stat_loess_bscore_zplate      =denv.k.hits.tstat.norm.loess.bsc.zplate$GeneSymbol,
                                  T_stat_loess_bscore_zctrl       =denv.k.hits.tstat.norm.loess.bsc.zctrl$GeneSymbol,
                                  T_stat_loess_bscore_zctrl_scram =denv.k.hits.tstat.norm.loess.bsc.zctrl.scram$GeneSymbol,
                                  H_stat_zplate                   =denv.k.hits.hstat.norm.zplate$GeneSymbol,
                                  H_stat_zctrl                    =denv.k.hits.hstat.norm.zctrl$GeneSymbol,
                                  H_stat_zctrl_scram              =denv.k.hits.hstat.norm.zctrl.scram$GeneSymbol,
                                  H_stat_bsc_zplate               =denv.k.hits.hstat.norm.bsc.zplate$GeneSymbol,
                                  H_stat_bsc_zctrl                =denv.k.hits.hstat.norm.bsc.zctrl$GeneSymbol,
                                  H_stat_bsc_zplate_scram         =denv.k.hits.hstat.norm.bsc.zctrl.scram$GeneSymbol,
                                  H_stat_loess_bscore_zplate      =denv.k.hits.hstat.norm.loess.bsc.zplate$GeneSymbol,
                                  H_stat_loess_bscore_zctrl       =denv.k.hits.hstat.norm.loess.bsc.zctrl$GeneSymbol,
                                  H_stat_loess_bscore_zctrl_scram =denv.k.hits.hstat.norm.loess.bsc.zctrl.scram$GeneSymbol ))

jac.kegg.hits <-  set.concordance(list(T_stat_zplate                   =kegg.mapping(denv.k.hits.tstat.norm.zplate$Entrez),
                                       T_stat_zctrl                    =kegg.mapping(denv.k.hits.tstat.norm.zctrl$Entrez),
                                       T_stat_zctrl_scram              =kegg.mapping(denv.k.hits.tstat.norm.zctrl.scram$Entrez),
                                       T_stat_bsc_zplate               =kegg.mapping(denv.k.hits.tstat.norm.bsc.zplate$Entrez),
                                       T_stat_bsc_zctrl                =kegg.mapping(denv.k.hits.tstat.norm.bsc.zctrl$Entrez),
                                       T_stat_bsc_zplate_scram         =kegg.mapping(denv.k.hits.tstat.norm.bsc.zctrl.scram$Entrez),
                                       T_stat_loess_bscore_zplate      =kegg.mapping(denv.k.hits.tstat.norm.loess.bsc.zplate$Entrez),
                                       T_stat_loess_bscore_zctrl       =kegg.mapping(denv.k.hits.tstat.norm.loess.bsc.zctrl$Entrez),
                                       T_stat_loess_bscore_zctrl_scram =kegg.mapping(denv.k.hits.tstat.norm.loess.bsc.zctrl.scram$Entrez),
                                       H_stat_zplate                   =kegg.mapping(denv.k.hits.hstat.norm.zplate$Entrez),
                                       H_stat_zctrl                    =kegg.mapping(denv.k.hits.hstat.norm.zctrl$Entrez),
                                       H_stat_zctrl_scram              =kegg.mapping(denv.k.hits.hstat.norm.zctrl.scram$Entrez),
                                       H_stat_bsc_zplate               =kegg.mapping(denv.k.hits.hstat.norm.bsc.zplate$Entrez),
                                       H_stat_bsc_zctrl                =kegg.mapping(denv.k.hits.hstat.norm.bsc.zctrl$Entrez),
                                       H_stat_bsc_zplate_scram         =kegg.mapping(denv.k.hits.hstat.norm.bsc.zctrl.scram$Entrez),
                                       H_stat_loess_bscore_zplate      =kegg.mapping(denv.k.hits.hstat.norm.loess.bsc.zplate$Entrez),
                                       H_stat_loess_bscore_zctrl       =kegg.mapping(denv.k.hits.hstat.norm.loess.bsc.zctrl$Entrez),
                                       H_stat_loess_bscore_zctrl_scram =kegg.mapping(denv.k.hits.hstat.norm.loess.bsc.zctrl.scram$Entrez)))

found.genes <- sort(table(c(denv.k.hits.tstat.norm.zplate$GeneSymbol,
                            denv.k.hits.tstat.norm.zctrl$GeneSymbol,
                            denv.k.hits.tstat.norm.zctrl.scram$GeneSymbol,
                            denv.k.hits.tstat.norm.bsc.zplate$GeneSymbol,
                            denv.k.hits.tstat.norm.bsc.zctrl$GeneSymbol,
                            denv.k.hits.tstat.norm.bsc.zctrl.scram$GeneSymbol,
                            denv.k.hits.tstat.norm.loess.bsc.zplate$GeneSymbol,
                            denv.k.hits.tstat.norm.loess.bsc.zctrl$GeneSymbol,
                            denv.k.hits.tstat.norm.loess.bsc.zctrl.scram$GeneSymbol,
                            denv.k.hits.hstat.norm.zplate$GeneSymbol,
                            denv.k.hits.hstat.norm.zctrl$GeneSymbol,
                            denv.k.hits.hstat.norm.zctrl.scram$GeneSymbol,
                            denv.k.hits.hstat.norm.bsc.zplate$GeneSymbol,
                            denv.k.hits.hstat.norm.bsc.zctrl$GeneSymbol,
                            denv.k.hits.hstat.norm.bsc.zctrl.scram$GeneSymbol,
                            denv.k.hits.hstat.norm.loess.bsc.zplate$GeneSymbol,
                            denv.k.hits.hstat.norm.loess.bsc.zctrl$GeneSymbol,
                            denv.k.hits.hstat.norm.loess.bsc.zctrl.scram$GeneSymbol)),
                    decreasing=T)
found.entrez <- sort(table(c(denv.k.hits.tstat.norm.zplate$Entrez,
                             denv.k.hits.tstat.norm.zctrl$Entrez,
                             denv.k.hits.tstat.norm.zctrl.scram$Entrez,
                             denv.k.hits.tstat.norm.bsc.zplate$Entrez,
                             denv.k.hits.tstat.norm.bsc.zctrl$Entrez,
                             denv.k.hits.tstat.norm.bsc.zctrl.scram$Entrez,
                             denv.k.hits.tstat.norm.loess.bsc.zplate$Entrez,
                             denv.k.hits.tstat.norm.loess.bsc.zctrl$Entrez,
                             denv.k.hits.tstat.norm.loess.bsc.zctrl.scram$Entrez,
                             denv.k.hits.hstat.norm.zplate$Entrez,
                             denv.k.hits.hstat.norm.zctrl$Entrez,
                             denv.k.hits.hstat.norm.zctrl.scram$Entrez,
                             denv.k.hits.hstat.norm.bsc.zplate$Entrez,
                             denv.k.hits.hstat.norm.bsc.zctrl$Entrez,
                             denv.k.hits.hstat.norm.bsc.zctrl.scram$Entrez,
                             denv.k.hits.hstat.norm.loess.bsc.zplate$Entrez,
                             denv.k.hits.hstat.norm.loess.bsc.zctrl$Entrez,
                             denv.k.hits.hstat.norm.loess.bsc.zctrl.scram$Entrez)),
                    decreasing=T)
# ORA
denv.k.ora.tstat.norm.zplate                 <- ora(denv.k.hits.tstat.norm.zplate$Entrez,                denv.k.tstat.norm.zplate$Entrez,                db="go")
denv.k.ora.tstat.norm.zctrl                  <- ora(denv.k.hits.tstat.norm.zctrl$Entrez,                 denv.k.tstat.norm.zctrl$Entrez,                 db="go")
denv.k.ora.tstat.norm.zctrl.scram            <- ora(denv.k.hits.tstat.norm.zctrl.scram$Entrez,           denv.k.tstat.norm.zctrl.scram$Entrez,           db="go")
denv.k.ora.tstat.norm.bsc.zplate             <- ora(denv.k.hits.tstat.norm.bsc.zplate$Entrez,            denv.k.tstat.norm.bsc.zplate$Entrez,            db="go")
denv.k.ora.tstat.norm.bsc.zctrl              <- ora(denv.k.hits.tstat.norm.bsc.zctrl$Entrez,             denv.k.tstat.norm.bsc.zctrl$Entrez,             db="go")
denv.k.ora.tstat.norm.bsc.zctrl.scram        <- ora(denv.k.hits.tstat.norm.bsc.zctrl.scram$Entrez,       denv.k.tstat.norm.bsc.zctrl.scram$Entrez,       db="go")
denv.k.ora.tstat.norm.loess.bsc.zplate       <- ora(denv.k.hits.tstat.norm.loess.bsc.zplate$Entrez,      denv.k.tstat.norm.loess.bsc.zplate$Entrez,      db="go")
denv.k.ora.tstat.norm.loess.bsc.zctrl        <- ora(denv.k.hits.tstat.norm.loess.bsc.zctrl$Entrez,       denv.k.tstat.norm.loess.bsc.zctrl$Entrez,       db="go")
denv.k.ora.tstat.norm.loess.bsc.zctrl.scram  <- ora(denv.k.hits.tstat.norm.loess.bsc.zctrl.scram$Entrez, denv.k.tstat.norm.loess.bsc.zctrl.scram$Entrez, db="go")
denv.k.ora.hstat.norm.zplate                 <- ora(denv.k.hits.hstat.norm.zplate$Entrez,                denv.k.hstat.norm.zplate$Entrez,                db="go")
denv.k.ora.hstat.norm.zctrl                  <- ora(denv.k.hits.hstat.norm.zctrl$Entrez,                 denv.k.hstat.norm.zctrl$Entrez,                 db="go")
denv.k.ora.hstat.norm.zctrl.scram            <- ora(denv.k.hits.hstat.norm.zctrl.scram$Entrez,           denv.k.hstat.norm.zctrl.scram$Entrez,           db="go")
denv.k.ora.hstat.norm.bsc.zplate             <- ora(denv.k.hits.hstat.norm.bsc.zplate$Entrez,            denv.k.hstat.norm.bsc.zplate$Entrez,            db="go")
denv.k.ora.hstat.norm.bsc.zctrl              <- ora(denv.k.hits.hstat.norm.bsc.zctrl$Entrez,             denv.k.hstat.norm.bsc.zctrl$Entrez,             db="go")
denv.k.ora.hstat.norm.bsc.zctrl.scram        <- ora(denv.k.hits.hstat.norm.bsc.zctrl.scram$Entrez,       denv.k.hstat.norm.bsc.zctrl.scram$Entrez,       db="go")
denv.k.ora.hstat.norm.loess.bsc.zplate       <- ora(denv.k.hits.hstat.norm.loess.bsc.zplate$Entrez,      denv.k.hstat.norm.loess.bsc.zplate$Entrez,      db="go")
denv.k.ora.hstat.norm.loess.bsc.zctrl        <- ora(denv.k.hits.hstat.norm.loess.bsc.zctrl$Entrez,       denv.k.hstat.norm.loess.bsc.zctrl$Entrez,       db="go")
denv.k.ora.hstat.norm.loess.bsc.zctrl.scram  <- ora(denv.k.hits.hstat.norm.loess.bsc.zctrl.scram$Entrez, denv.k.hstat.norm.loess.bsc.zctrl.scram$Entrez, db="go")

all.genes.go.ora   <- ora(unique(as.integer(names(found.entrez))), denv.k.summ.loess.bsc.zctrl.scram$Entrez, db="go")
all.genes.kegg.ora <- ora(unique(as.integer(names(found.entrez))), denv.k.summ.loess.bsc.zctrl.scram$Entrez, db="kegg")

jac.ora <- set.concordance(list( T_stat_zplate                    =denv.k.ora.tstat.norm.zplate$summary$Term,
                                 T_stat_zctrl                     =denv.k.ora.tstat.norm.zctrl$summary$Term,
                                 T_stat_zctrl_scram               =denv.k.ora.tstat.norm.zctrl.scram$summary$Term,
                                 T_stat_bsc_zplate                =denv.k.ora.tstat.norm.bsc.zplate$summary$Term,
                                 T_stat_bsc_zctrl                 =denv.k.ora.tstat.norm.bsc.zctrl$summary$Term,
                                 T_stat_bsc_zctrl_scram           =denv.k.ora.tstat.norm.bsc.zctrl.scram$summary$Term,
                                 T_stat_loess_bsc_zplate          =denv.k.ora.tstat.norm.loess.bsc.zplate$summary$Term,
                                 T_stat_loess_bsc_zctrl           =denv.k.ora.tstat.norm.loess.bsc.zctrl$summary$Term,
                                 T_stat_loess_bsc_zctrl_scram     =denv.k.ora.tstat.norm.loess.bsc.zctrl.scram$summary$Term,
                                 H_stat_zplate                    =denv.k.ora.hstat.norm.zplate$summary$Term,
                                 H_stat_zctrl                     =denv.k.ora.hstat.norm.zctrl$summary$Term,
                                 H_stat_zctrl_scram               =denv.k.ora.hstat.norm.zctrl.scram$summary$Term,
                                 H_stat_bsc_zplate                =denv.k.ora.hstat.norm.bsc.zplate$summary$Term,
                                 H_stat_bsc_zctrl                 =denv.k.ora.hstat.norm.bsc.zctrl$summary$Term,
                                 H_stat_bsc_zctrl_scram           =denv.k.ora.hstat.norm.bsc.zctrl.scram$summary$Term,
                                 H_stat_loess_bsc_zplate          =denv.k.ora.hstat.norm.loess.bsc.zplate$summary$Term,
                                 H_stat_loess_bsc_zctrl           =denv.k.ora.hstat.norm.loess.bsc.zctrl$summary$Term,
                                 H_stat_loess_bsc_zctrl_scram     =denv.k.ora.hstat.norm.loess.bsc.zctrl.scram$summary$Term))

used.go.terms <- sort(table(c( denv.k.ora.tstat.norm.zplate$summary$Term,
                               denv.k.ora.tstat.norm.zctrl$summary$Term,
                               denv.k.ora.tstat.norm.zctrl.scram$summary$Term,
                               denv.k.ora.tstat.norm.bsc.zplate$summary$Term,
                               denv.k.ora.tstat.norm.bsc.zctrl$summary$Term,
                               denv.k.ora.tstat.norm.bsc.zctrl.scram$summary$Term,
                               denv.k.ora.tstat.norm.loess.bsc.zplate$summary$Term,
                               denv.k.ora.tstat.norm.loess.bsc.zctrl$summary$Term,
                               denv.k.ora.tstat.norm.loess.bsc.zctrl.scram$summary$Term,
                               denv.k.ora.hstat.norm.zplate$summary$Term,
                               denv.k.ora.hstat.norm.zctrl$summary$Term,
                               denv.k.ora.hstat.norm.zctrl.scram$summary$Term,
                               denv.k.ora.hstat.norm.bsc.zplate$summary$Term,
                               denv.k.ora.hstat.norm.bsc.zctrl$summary$Term,
                               denv.k.ora.hstat.norm.bsc.zctrl.scram$summary$Term,
                               denv.k.ora.hstat.norm.loess.bsc.zplate$summary$Term,
                               denv.k.ora.hstat.norm.loess.bsc.zctrl$summary$Term,
                               denv.k.ora.hstat.norm.loess.bsc.zctrl.scram$summary$Term)),
                   decreasing=T)

# write all to file
df.used.terms <- data.frame(Number=as.vector(used.go.terms), Name=names(used.terms))
df.found.genes <- data.frame(Number=as.vector(found.genes), Name=names(found.genes))
hit.p1 <- plot(jac.hits)
ggsave(hit.p1, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/single_virus_analysis/denv_kinome_analysis_jaccard_hits.pdf")
ora.p1 <- plot(jac.ora)
ggsave(ora.p1, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/single_virus_analysis/denv_kinome_analysis_jaccard_ora.pdf")
kegg.p1 <- plot(jac.kegg.hits)
ggsave(kegg.p1, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/single_virus_analysis/denv_kinome_analysis_jaccard_kegg.pdf")

write.table(df.used.terms,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/denv_kinome_analysis_table_ora.tsv",
            quote=F, row.names=F, sep="\t")

write.table(df.found.genes,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/denv_kinome_analysis_table_genes.tsv",
            quote=F, row.names=F, sep="\t")

write.table(all.genes.go.ora$summary,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/denv_kinome_analysis_table_combined_gene_list_go_ora.tsv",
            quote=F, row.names=F, sep="\t")


write.table(all.genes.kegg.ora$summary,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/denv_kinome_analysis_table_combined_gene_list_kegg_ora.tsv",
            quote=F, row.names=F, sep="\t")
