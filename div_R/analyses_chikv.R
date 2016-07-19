library(data.table, quietly=T)
library(dplyr, quietly=T)
library(tidyr, quietly=T)
library(reshape2, quietly=T)
library(svd)

rm(list=setdiff(ls(), "rnai.screen.raw"))
load("inst/extdata/rnai.screen.raw.rda")

chikv.raw        <- filter(rnai.screen.raw, Virus=="CHIKV")
chikv.quality.raw <- quality(chikv.raw)
plot(chikv.quality.raw)

chikv.norm       <- svd::preprocess(chikv.raw, normalize=c("log", "background", "robust-z.score"),
                                    method="mean",
                                    background.column=12,
                                    z.score.level="plate",
                                    poc.ctrl="Scrambled",
                                    drop=T,
                                    rm.cytotoxic="Scrambled")
chikv.hyper.ana <- hyperstatistic(chikv.norm)
chikv.priorit <- prioritize(chikv.hyper.ana, hit.ratio=.66, readout.threshold=2)


hist(chikv.raw)
hist(chikv.norm)

chikv.plates.raw <- plates(chikv.raw)
chikv.plates.norm <- plates(chikv.norm)
sam <- sample(length(chikv.plates.raw), 1)
svd:::.multiplot(
  plot(chikv.plates.raw[[sam]], show.gene.names=T),
  plot(chikv.plates.norm[[sam]], show.gene.names=T))

chikv.replicates.raw <- replicates(chikv.raw)
chikv.replicates.norm <- replicates(chikv.norm)
plot(chikv.replicates.raw)
plot(chikv.replicates.norm)

plot(chikv.hyper.ana, readout.thresh=2, sig.thresh=.2)
plot(chikv.priorit)

