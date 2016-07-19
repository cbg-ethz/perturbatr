library(data.table, quietly=T)
library(dplyr, quietly=T)
library(tidyr, quietly=T)
library(reshape2, quietly=T)

rm(list = setdiff(ls(), lsf.str()))

load("inst/extdata/rnai.screen.raw.rda")


### hits for chikv using tstatistics and hyperstatistics
# get hits
hcv.k.dat              <- filter(rnai.screen.raw, Virus=="HCV", Screen=="Kinome")
# normalize
hcv.k.norm.zplate                <- normalize(hcv.k.dat, normalize=c("robust-z.score"), level="plate")
hcv.k.norm.zctrl                 <- normalize(hcv.k.dat, normalize=c("robust-z.score"), level="control")
hcv.k.norm.zctrl.scram           <- normalize(hcv.k.dat, normalize=c("robust-z.score"), level="control", ctrl="Scrambled")
hcv.k.norm.bsc.zplate            <- normalize(hcv.k.dat, normalize=c("b.score", "robust-z.score"), level="plate")
hcv.k.norm.bsc.zctrl             <- normalize(hcv.k.dat, normalize=c("b.score", "robust-z.score"), level="control")
hcv.k.norm.bsc.zctrl.scram       <- normalize(hcv.k.dat, normalize=c("b.score", "robust-z.score"), level="control", ctrl="Scrambled")
hcv.k.norm.loess.bsc.zplate      <- normalize(hcv.k.dat, normalize=c("loess", "b.score", "robust-z.score"), level="plate")
hcv.k.norm.loess.bsc.zctrl       <- normalize(hcv.k.dat, normalize=c("loess", "b.score", "robust-z.score"), level="control")
hcv.k.norm.loess.bsc.zctrl.scram <- normalize(hcv.k.dat, normalize=c("loess", "b.score", "robust-z.score"), level="control", ctrl="Scrambled")
# summarize
hcv.k.summ.zplate                <- summarize(hcv.k.norm.zplate)
hcv.k.summ.zctrl                 <- summarize(hcv.k.norm.zctrl)
hcv.k.summ.zctrl.scram           <- summarize(hcv.k.norm.zctrl.scram)
hcv.k.summ.bsc.zplate            <- summarize(hcv.k.norm.bsc.zplate)
hcv.k.summ.bsc.zctrl             <- summarize(hcv.k.norm.bsc.zctrl)
hcv.k.summ.bsc.zctrl.scram       <- summarize(hcv.k.norm.bsc.zctrl.scram)
hcv.k.summ.loess.bsc.zplate      <- summarize(hcv.k.norm.loess.bsc.zplate)
hcv.k.summ.loess.bsc.zctrl       <- summarize(hcv.k.norm.loess.bsc.zctrl)
hcv.k.summ.loess.bsc.zctrl.scram <- summarize(hcv.k.norm.loess.bsc.zctrl.scram)
# analysis
hcv.k.tstat.norm.zplate                   <- tstatistic(hcv.k.summ.zplate,                padjust="BH")
hcv.k.tstat.norm.zctrl                    <- tstatistic(hcv.k.summ.zctrl,                 padjust="BH", mu="control")
hcv.k.tstat.norm.zctrl.scram              <- tstatistic(hcv.k.summ.zctrl.scram,           padjust="BH", mu="control", ctrl="Scrambled")
hcv.k.tstat.norm.bsc.zplate               <- tstatistic(hcv.k.summ.bsc.zplate,            padjust="BH")
hcv.k.tstat.norm.bsc.zctrl                <- tstatistic(hcv.k.summ.bsc.zctrl,             padjust="BH", mu="control")
hcv.k.tstat.norm.bsc.zctrl.scram          <- tstatistic(hcv.k.summ.bsc.zctrl.scram,       padjust="BH", mu="control", ctrl="Scrambled")
hcv.k.tstat.norm.loess.bsc.zplate         <- tstatistic(hcv.k.summ.loess.bsc.zplate,      padjust="BH")
hcv.k.tstat.norm.loess.bsc.zctrl          <- tstatistic(hcv.k.summ.loess.bsc.zctrl,       padjust="BH", mu="control")
hcv.k.tstat.norm.loess.bsc.zctrl.scram    <- tstatistic(hcv.k.summ.loess.bsc.zctrl.scram, padjust="BH", mu="control", ctrl="Scrambled")
hcv.k.hstat.norm.zplate                   <- hyperstatistic(hcv.k.summ.zplate)
hcv.k.hstat.norm.zctrl                    <- hyperstatistic(hcv.k.summ.zctrl)
hcv.k.hstat.norm.zctrl.scram              <- hyperstatistic(hcv.k.summ.zctrl.scram)
hcv.k.hstat.norm.bsc.zplate               <- hyperstatistic(hcv.k.summ.bsc.zplate)
hcv.k.hstat.norm.bsc.zctrl                <- hyperstatistic(hcv.k.summ.bsc.zctrl)
hcv.k.hstat.norm.bsc.zctrl.scram          <- hyperstatistic(hcv.k.summ.bsc.zctrl.scram)
hcv.k.hstat.norm.loess.bsc.zplate         <- hyperstatistic(hcv.k.summ.loess.bsc.zplate)
hcv.k.hstat.norm.loess.bsc.zctrl          <- hyperstatistic(hcv.k.summ.loess.bsc.zctrl)
hcv.k.hstat.norm.loess.bsc.zctrl.scram    <- hyperstatistic(hcv.k.summ.loess.bsc.zctrl.scram)
# hit calling
hcv.k.hits.tstat.norm.zplate                 <- prioritize(hcv.k.tstat.norm.zplate,                hit.ratio=.6)
hcv.k.hits.tstat.norm.zctrl                  <- prioritize(hcv.k.tstat.norm.zctrl,                 hit.ratio=.6)
hcv.k.hits.tstat.norm.zctrl.scram            <- prioritize(hcv.k.tstat.norm.zctrl.scram,           hit.ratio=.6)
hcv.k.hits.tstat.norm.bsc.zplate             <- prioritize(hcv.k.tstat.norm.bsc.zplate,            hit.ratio=.6)
hcv.k.hits.tstat.norm.bsc.zctrl              <- prioritize(hcv.k.tstat.norm.bsc.zctrl,             hit.ratio=.6)
hcv.k.hits.tstat.norm.bsc.zctrl.scram        <- prioritize(hcv.k.tstat.norm.bsc.zctrl.scram,       hit.ratio=.6)
hcv.k.hits.tstat.norm.loess.bsc.zplate       <- prioritize(hcv.k.tstat.norm.loess.bsc.zplate,      hit.ratio=.6)
hcv.k.hits.tstat.norm.loess.bsc.zctrl        <- prioritize(hcv.k.tstat.norm.loess.bsc.zctrl,       hit.ratio=.6)
hcv.k.hits.tstat.norm.loess.bsc.zctrl.scram  <- prioritize(hcv.k.tstat.norm.loess.bsc.zctrl.scram, hit.ratio=.6, readout.thresh=2)
hcv.k.hits.hstat.norm.zplate                 <- prioritize(hcv.k.hstat.norm.zplate,                hit.ratio=.6)
hcv.k.hits.hstat.norm.zctrl                  <- prioritize(hcv.k.hstat.norm.zctrl,                 hit.ratio=.6)
hcv.k.hits.hstat.norm.zctrl.scram            <- prioritize(hcv.k.hstat.norm.zctrl.scram,           hit.ratio=.6)
hcv.k.hits.hstat.norm.bsc.zplate             <- prioritize(hcv.k.hstat.norm.bsc.zplate,            hit.ratio=.6)
hcv.k.hits.hstat.norm.bsc.zctrl              <- prioritize(hcv.k.hstat.norm.bsc.zctrl,             hit.ratio=.6)
hcv.k.hits.hstat.norm.bsc.zctrl.scram        <- prioritize(hcv.k.hstat.norm.bsc.zctrl.scram,       hit.ratio=.6)
hcv.k.hits.hstat.norm.loess.bsc.zplate       <- prioritize(hcv.k.hstat.norm.loess.bsc.zplate,      hit.ratio=.6)
hcv.k.hits.hstat.norm.loess.bsc.zctrl        <- prioritize(hcv.k.hstat.norm.loess.bsc.zctrl,       hit.ratio=.6)
hcv.k.hits.hstat.norm.loess.bsc.zctrl.scram  <- prioritize(hcv.k.hstat.norm.loess.bsc.zctrl.scram, hit.ratio=.6, readout.thresh=2)

# get number of hits
jac.hits <- set.concordance(list( T_stat_zplate                   =hcv.k.hits.tstat.norm.zplate$GeneSymbol,
                                  T_stat_zctrl                    =hcv.k.hits.tstat.norm.zctrl$GeneSymbol,
                                  T_stat_zctrl_scram              =hcv.k.hits.tstat.norm.zctrl.scram$GeneSymbol,
                                  T_stat_bsc_zplate               =hcv.k.hits.tstat.norm.bsc.zplate$GeneSymbol,
                                  T_stat_bsc_zctrl                =hcv.k.hits.tstat.norm.bsc.zctrl$GeneSymbol,
                                  T_stat_bsc_zplate_scram         =hcv.k.hits.tstat.norm.bsc.zctrl.scram$GeneSymbol,
                                  T_stat_loess_bscore_zplate      =hcv.k.hits.tstat.norm.loess.bsc.zplate$GeneSymbol,
                                  T_stat_loess_bscore_zctrl       =hcv.k.hits.tstat.norm.loess.bsc.zctrl$GeneSymbol,
                                  T_stat_loess_bscore_zctrl_scram =hcv.k.hits.tstat.norm.loess.bsc.zctrl.scram$GeneSymbol,
                                  H_stat_zplate                   =hcv.k.hits.hstat.norm.zplate$GeneSymbol,
                                  H_stat_zctrl                    =hcv.k.hits.hstat.norm.zctrl$GeneSymbol,
                                  H_stat_zctrl_scram              =hcv.k.hits.hstat.norm.zctrl.scram$GeneSymbol,
                                  H_stat_bsc_zplate               =hcv.k.hits.hstat.norm.bsc.zplate$GeneSymbol,
                                  H_stat_bsc_zctrl                =hcv.k.hits.hstat.norm.bsc.zctrl$GeneSymbol,
                                  H_stat_bsc_zplate_scram         =hcv.k.hits.hstat.norm.bsc.zctrl.scram$GeneSymbol,
                                  H_stat_loess_bscore_zplate      =hcv.k.hits.hstat.norm.loess.bsc.zplate$GeneSymbol,
                                  H_stat_loess_bscore_zctrl       =hcv.k.hits.hstat.norm.loess.bsc.zctrl$GeneSymbol,
                                  H_stat_loess_bscore_zctrl_scram =hcv.k.hits.hstat.norm.loess.bsc.zctrl.scram$GeneSymbol ))

jac.kegg.hits <- set.concordance(list( T_stat_zplate              =kegg.mapping(hcv.k.hits.tstat.norm.zplate$Entrez),
                                  T_stat_zctrl                    =kegg.mapping(hcv.k.hits.tstat.norm.zctrl$Entrez),
                                  T_stat_zctrl_scram              =kegg.mapping(hcv.k.hits.tstat.norm.zctrl.scram$Entrez),
                                  T_stat_bsc_zplate               =kegg.mapping(hcv.k.hits.tstat.norm.bsc.zplate$Entrez),
                                  T_stat_bsc_zctrl                =kegg.mapping(hcv.k.hits.tstat.norm.bsc.zctrl$Entrez),
                                  T_stat_bsc_zplate_scram         =kegg.mapping(hcv.k.hits.tstat.norm.bsc.zctrl.scram$Entrez),
                                  T_stat_loess_bscore_zplate      =kegg.mapping(hcv.k.hits.tstat.norm.loess.bsc.zplate$Entrez),
                                  T_stat_loess_bscore_zctrl       =kegg.mapping(hcv.k.hits.tstat.norm.loess.bsc.zctrl$Entrez),
                                  T_stat_loess_bscore_zctrl_scram =kegg.mapping(hcv.k.hits.tstat.norm.loess.bsc.zctrl.scram$Entrez),
                                  H_stat_zplate                   =kegg.mapping(hcv.k.hits.hstat.norm.zplate$Entrez),
                                  H_stat_zctrl                    =kegg.mapping(hcv.k.hits.hstat.norm.zctrl$Entrez),
                                  H_stat_zctrl_scram              =kegg.mapping(hcv.k.hits.hstat.norm.zctrl.scram$Entrez),
                                  H_stat_bsc_zplate               =kegg.mapping(hcv.k.hits.hstat.norm.bsc.zplate$Entrez),
                                  H_stat_bsc_zctrl                =kegg.mapping(hcv.k.hits.hstat.norm.bsc.zctrl$Entrez),
                                  H_stat_bsc_zplate_scram         =kegg.mapping(hcv.k.hits.hstat.norm.bsc.zctrl.scram$Entrez),
                                  H_stat_loess_bscore_zplate      =kegg.mapping(hcv.k.hits.hstat.norm.loess.bsc.zplate$Entrez),
                                  H_stat_loess_bscore_zctrl       =kegg.mapping(hcv.k.hits.hstat.norm.loess.bsc.zctrl$Entrez),
                                  H_stat_loess_bscore_zctrl_scram =kegg.mapping(hcv.k.hits.hstat.norm.loess.bsc.zctrl.scram$Entrez) ))

found.genes <- sort(table(c(hcv.k.hits.tstat.norm.zplate$GeneSymbol,
                            hcv.k.hits.tstat.norm.zctrl$GeneSymbol,
                            hcv.k.hits.tstat.norm.zctrl.scram$GeneSymbol,
                            hcv.k.hits.tstat.norm.bsc.zplate$GeneSymbol,
                            hcv.k.hits.tstat.norm.bsc.zctrl$GeneSymbol,
                            hcv.k.hits.tstat.norm.bsc.zctrl.scram$GeneSymbol,
                            hcv.k.hits.tstat.norm.loess.bsc.zplate$GeneSymbol,
                            hcv.k.hits.tstat.norm.loess.bsc.zctrl$GeneSymbol,
                            hcv.k.hits.tstat.norm.loess.bsc.zctrl.scram$GeneSymbol,
                            hcv.k.hits.hstat.norm.zplate$GeneSymbol,
                            hcv.k.hits.hstat.norm.zctrl$GeneSymbol,
                            hcv.k.hits.hstat.norm.zctrl.scram$GeneSymbol,
                            hcv.k.hits.hstat.norm.bsc.zplate$GeneSymbol,
                            hcv.k.hits.hstat.norm.bsc.zctrl$GeneSymbol,
                            hcv.k.hits.hstat.norm.bsc.zctrl.scram$GeneSymbol,
                            hcv.k.hits.hstat.norm.loess.bsc.zplate$GeneSymbol,
                            hcv.k.hits.hstat.norm.loess.bsc.zctrl$GeneSymbol,
                            hcv.k.hits.hstat.norm.loess.bsc.zctrl.scram$GeneSymbol)),
                    decreasing=T)

found.entrez <- sort(table(c(hcv.k.hits.tstat.norm.zplate$Entrez,
                             hcv.k.hits.tstat.norm.zctrl$Entrez,
                             hcv.k.hits.tstat.norm.zctrl.scram$Entrez,
                             hcv.k.hits.tstat.norm.bsc.zplate$Entrez,
                             hcv.k.hits.tstat.norm.bsc.zctrl$Entrez,
                             hcv.k.hits.tstat.norm.bsc.zctrl.scram$Entrez,
                             hcv.k.hits.tstat.norm.loess.bsc.zplate$Entrez,
                             hcv.k.hits.tstat.norm.loess.bsc.zctrl$Entrez,
                             hcv.k.hits.tstat.norm.loess.bsc.zctrl.scram$Entrez,
                             hcv.k.hits.hstat.norm.zplate$Entrez,
                             hcv.k.hits.hstat.norm.zctrl$Entrez,
                             hcv.k.hits.hstat.norm.zctrl.scram$Entrez,
                             hcv.k.hits.hstat.norm.bsc.zplate$Entrez,
                             hcv.k.hits.hstat.norm.bsc.zctrl$Entrez,
                             hcv.k.hits.hstat.norm.bsc.zctrl.scram$Entrez,
                             hcv.k.hits.hstat.norm.loess.bsc.zplate$Entrez,
                             hcv.k.hits.hstat.norm.loess.bsc.zctrl$Entrez,
                             hcv.k.hits.hstat.norm.loess.bsc.zctrl.scram$Entrez)),
                     decreasing=T)

# ORA
hcv.k.ora.tstat.norm.zplate                 <- ora(hcv.k.hits.tstat.norm.zplate$Entrez,                hcv.k.tstat.norm.zplate$Entrez,                db="go")
hcv.k.ora.tstat.norm.zctrl                  <- ora(hcv.k.hits.tstat.norm.zctrl$Entrez,                 hcv.k.tstat.norm.zctrl$Entrez,                 db="go")
hcv.k.ora.tstat.norm.zctrl.scram            <- ora(hcv.k.hits.tstat.norm.zctrl.scram$Entrez,           hcv.k.tstat.norm.zctrl.scram$Entrez,           db="go")
hcv.k.ora.tstat.norm.bsc.zplate             <- ora(hcv.k.hits.tstat.norm.bsc.zplate$Entrez,            hcv.k.tstat.norm.bsc.zplate$Entrez,            db="go")
hcv.k.ora.tstat.norm.bsc.zctrl              <- ora(hcv.k.hits.tstat.norm.bsc.zctrl$Entrez,             hcv.k.tstat.norm.bsc.zctrl$Entrez,             db="go")
hcv.k.ora.tstat.norm.bsc.zctrl.scram        <- ora(hcv.k.hits.tstat.norm.bsc.zctrl.scram$Entrez,       hcv.k.tstat.norm.bsc.zctrl.scram$Entrez,       db="go")
hcv.k.ora.tstat.norm.loess.bsc.zplate       <- ora(hcv.k.hits.tstat.norm.loess.bsc.zplate$Entrez,      hcv.k.tstat.norm.loess.bsc.zplate$Entrez,      db="go")
hcv.k.ora.tstat.norm.loess.bsc.zctrl        <- ora(hcv.k.hits.tstat.norm.loess.bsc.zctrl$Entrez,       hcv.k.tstat.norm.loess.bsc.zctrl$Entrez,       db="go")
hcv.k.ora.tstat.norm.loess.bsc.zctrl.scram  <- ora(hcv.k.hits.tstat.norm.loess.bsc.zctrl.scram$Entrez, hcv.k.tstat.norm.loess.bsc.zctrl.scram$Entrez, db="go")
hcv.k.ora.hstat.norm.zplate                 <- ora(hcv.k.hits.hstat.norm.zplate$Entrez,                hcv.k.hstat.norm.zplate$Entrez,                db="go")
hcv.k.ora.hstat.norm.zctrl                  <- ora(hcv.k.hits.hstat.norm.zctrl$Entrez,                 hcv.k.hstat.norm.zctrl$Entrez,                 db="go")
hcv.k.ora.hstat.norm.zctrl.scram            <- ora(hcv.k.hits.hstat.norm.zctrl.scram$Entrez,           hcv.k.hstat.norm.zctrl.scram$Entrez,           db="go")
hcv.k.ora.hstat.norm.bsc.zplate             <- ora(hcv.k.hits.hstat.norm.bsc.zplate$Entrez,            hcv.k.hstat.norm.bsc.zplate$Entrez,            db="go")
hcv.k.ora.hstat.norm.bsc.zctrl              <- ora(hcv.k.hits.hstat.norm.bsc.zctrl$Entrez,             hcv.k.hstat.norm.bsc.zctrl$Entrez,             db="go")
hcv.k.ora.hstat.norm.bsc.zctrl.scram        <- ora(hcv.k.hits.hstat.norm.bsc.zctrl.scram$Entrez,       hcv.k.hstat.norm.bsc.zctrl.scram$Entrez,       db="go")
hcv.k.ora.hstat.norm.loess.bsc.zplate       <- ora(hcv.k.hits.hstat.norm.loess.bsc.zplate$Entrez,      hcv.k.hstat.norm.loess.bsc.zplate$Entrez,      db="go")
hcv.k.ora.hstat.norm.loess.bsc.zctrl        <- ora(hcv.k.hits.hstat.norm.loess.bsc.zctrl$Entrez,       hcv.k.hstat.norm.loess.bsc.zctrl$Entrez,       db="go")
hcv.k.ora.hstat.norm.loess.bsc.zctrl.scram  <- ora(hcv.k.hits.hstat.norm.loess.bsc.zctrl.scram$Entrez, hcv.k.hstat.norm.loess.bsc.zctrl.scram$Entrez, db="go")

all.genes.go.ora   <- ora(unique(as.integer(names(found.entrez))), hcv.k.hstat.norm.bsc.zctrl.scram$Entrez, db="go")
all.genes.kegg.ora <- ora(unique(as.integer(names(found.entrez))), hcv.k.hstat.norm.bsc.zctrl.scram$Entrez, db="kegg")

jac.ora <- set.concordance(list( T_stat_zplate                    =hcv.k.ora.tstat.norm.zplate$summary$Term,
                                 T_stat_zctrl                     =hcv.k.ora.tstat.norm.zctrl$summary$Term,
                                 T_stat_zctrl_scram               =hcv.k.ora.tstat.norm.zctrl.scram$summary$Term,
                                 T_stat_bsc_zplate                =hcv.k.ora.tstat.norm.bsc.zplate$summary$Term,
                                 T_stat_bsc_zctrl                 =hcv.k.ora.tstat.norm.bsc.zctrl$summary$Term,
                                 T_stat_bsc_zctrl_scram           =hcv.k.ora.tstat.norm.bsc.zctrl.scram$summary$Term,
                                 T_stat_loess_bsc_zplate          =hcv.k.ora.tstat.norm.loess.bsc.zplate$summary$Term,
                                 T_stat_loess_bsc_zctrl           =hcv.k.ora.tstat.norm.loess.bsc.zctrl$summary$Term,
                                 T_stat_loess_bsc_zctrl_scram     =hcv.k.ora.tstat.norm.loess.bsc.zctrl.scram$summary$Term,
                                 H_stat_zplate                    =hcv.k.ora.hstat.norm.zplate$summary$Term,
                                 H_stat_zctrl                     =hcv.k.ora.hstat.norm.zctrl$summary$Term,
                                 H_stat_zctrl_scram               =hcv.k.ora.hstat.norm.zctrl.scram$summary$Term,
                                 H_stat_bsc_zplate                =hcv.k.ora.hstat.norm.bsc.zplate$summary$Term,
                                 H_stat_bsc_zctrl                 =hcv.k.ora.hstat.norm.bsc.zctrl$summary$Term,
                                 H_stat_bsc_zctrl_scram           =hcv.k.ora.hstat.norm.bsc.zctrl.scram$summary$Term,
                                 H_stat_loess_bsc_zplate          =hcv.k.ora.hstat.norm.loess.bsc.zplate$summary$Term,
                                 H_stat_loess_bsc_zctrl           =hcv.k.ora.hstat.norm.loess.bsc.zctrl$summary$Term,
                                 H_stat_loess_bsc_zctrl_scram     =hcv.k.ora.hstat.norm.loess.bsc.zctrl.scram$summary$Term))

used.terms <- sort(table(c( hcv.k.ora.tstat.norm.zplate$summary$Term,
                            hcv.k.ora.tstat.norm.zctrl$summary$Term,
                            hcv.k.ora.tstat.norm.zctrl.scram$summary$Term,
                            hcv.k.ora.tstat.norm.bsc.zplate$summary$Term,
                            hcv.k.ora.tstat.norm.bsc.zctrl$summary$Term,
                            hcv.k.ora.tstat.norm.bsc.zctrl.scram$summary$Term,
                            hcv.k.ora.tstat.norm.loess.bsc.zplate$summary$Term,
                            hcv.k.ora.tstat.norm.loess.bsc.zctrl$summary$Term,
                            hcv.k.ora.tstat.norm.loess.bsc.zctrl.scram$summary$Term,
                            hcv.k.ora.hstat.norm.zplate$summary$Term,
                            hcv.k.ora.hstat.norm.zctrl$summary$Term,
                            hcv.k.ora.hstat.norm.zctrl.scram$summary$Term,
                            hcv.k.ora.hstat.norm.bsc.zplate$summary$Term,
                            hcv.k.ora.hstat.norm.bsc.zctrl$summary$Term,
                            hcv.k.ora.hstat.norm.bsc.zctrl.scram$summary$Term,
                            hcv.k.ora.hstat.norm.loess.bsc.zplate$summary$Term,
                            hcv.k.ora.hstat.norm.loess.bsc.zctrl$summary$Term,
                            hcv.k.ora.hstat.norm.loess.bsc.zctrl.scram$summary$Term)),
                    decreasing=T)

# write all to file
df.used.terms <- data.frame(Number=as.vector(used.terms), Name=names(used.terms))
df.found.genes <- data.frame(Number=as.vector(found.genes), Name=names(found.genes))
hit.p1 <- plot(jac.hits)
ggsave(hit.p1, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/single_virus_analysis/hcv_kinome_analysis_jaccard_hits.pdf")
ora.p1 <- plot(jac.ora)
ggsave(ora.p1, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/single_virus_analysis/hcv_kinome_analysis_jaccard_ora.pdf")
kegg.p1 <- plot(jac.kegg.hits)
ggsave(ora.p1, filename="/Users/simondi/PHD/data/data/sysvirdrug/target_results/plots/single_virus_analysis/hcv_kinome_analysis_jaccard_kegg.pdf")

write.table(df.used.terms,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/hcv_kinome_analysis_table_ora.tsv",
            quote=F, row.names=F, sep="\t")

write.table(df.found.genes,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/hcv_kinome_analysis_table_genes.tsv",
            quote=F, row.names=F, sep="\t")

write.table(all.genes.go.ora$summary,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/hcv_kinome_analysis_table_combined_gene_list_go_ora.tsv",
            quote=F, row.names=F, sep="\t")

write.table(all.genes.kegg.ora$summary,
            file="/Users/simondi/PHD/data/data/sysvirdrug/target_results/tables/single_virus_analysis/hcv_kinome_analysis_table_combined_gene_list_kegg_ora.tsv",
            quote=F, row.names=F, sep="\t")
