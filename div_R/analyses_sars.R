library(data.table, quietly=T)
library(dplyr, quietly=T)
library(tidyr, quietly=T)
library(reshape2, quietly=T)
library(svd)

rm(list=setdiff(ls(), "rnai.screen.raw"))
load("inst/extdata/rnai.screen.raw.rda")

sars.raw        <- filter(rnai.screen.raw, Virus=="SARS")
sars.quality.raw <- quality(sars.raw)
plot(sars.quality.raw)

sars.norm       <- svd::preprocess(sars.raw, normalize=c("log", "background", "robust-z.score"),
                                    method="mean",
                                    background.column=12,
                                    z.score.level="plate",
                                    poc.ctrl="Scrambled",
                                    drop=T,
                                    rm.cytotoxic="Scrambled")
sars.hyper.ana <- hyperstatistic(sars.norm, level="gene")

sars.priorit <- prioritize(sars.hyper.ana, hit.ratio=.66, readout.threshold=2)


hist(sars.raw)
hist(sars.norm)

sars.plates.raw <- plates(sars.raw)
plot(sars.plates.raw, count=10)
sars.plates.norm <- plates(sars.norm)
sam <- sample(length(sars.plates.raw), 1)
svd:::.multiplot(
  plot(sars.plates.raw[[sam]], show.gene.names=T),
  plot(sars.plates.norm[[sam]], show.gene.names=T))

sars.replicates.raw <- replicates(sars.raw)
sars.replicates.norm <- replicates(sars.norm)
plot(sars.replicates.raw)
plot(sars.replicates.norm)

plot(sars.hyper.ana, readout.thresh=2, sig.thresh=.2)
plot(sars.priorit)

