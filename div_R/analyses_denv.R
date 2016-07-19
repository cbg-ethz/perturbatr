library(data.table, quietly=T)
library(dplyr, quietly=T)
library(dtplyr, quietly=T)
library(tidyr, quietly=T)
library(reshape2, quietly=T)
library(svd)

load("inst/extdata/rnai.screen.raw.rda")
rm(list=setdiff(ls(), "rnai.screen.raw"))

denv.raw        <- filter(rnai.screen.raw, Virus=="DENV", Screen=="Kinome")

denv.replicates.raw <- replicates(denv.raw)
plot(denv.replicates.raw)

denv.plates.raw <- plates(denv.raw)
plot(denv.plates.raw, count=10, show.gene.names=T)

denv.quality.raw <- quality(denv.raw)
plot(denv.quality.raw)

denv.norm       <- svd::preprocess(denv.raw, normalize=c("loess", "b.score",   "robust-z.score"),
                                  method="mean",
                                  z.score.level="plate",
                                  drop=F,
                                  rm.cytotoxic="Scrambled",
                                  #rm.outlier.wells="q",
                                  outlier.well.range=c(.1, .9))

denv.quality.norm <- quality(denv.norm)
plot(denv.quality.norm)

denv.replicates.norm <- replicates(denv.norm)
svd::cor(denv.replicates.norm) %>% data.table::melt() %>% summarize(m=mean(value))


denv.hyper.ana <- hyperstatistic(denv.norm, level="sirna")

denv.t.ana     <- tstatistic(denv.norm)

denv.priorit.h <- prioritize(denv.hyper.ana, hit.ratio=.5)


hist(denv.raw)
hist(denv.norm)

denv.plates.raw <- plates(denv.raw)
denv.plates.norm <- plates(denv.norm)
sam <- sample(length(denv.plates.raw), 1)
plot(denv.plates.raw, count=3,  show.controls=T)
plot(denv.plates.norm , count=3, show.controls=T)
svd:::.multiplot(
  plot(denv.plates.raw[[sam]], show.gene.names=T),
  plot(denv.plates.norm[[sam]], show.gene.names=T))

denv.replicates.raw <- replicates(denv.raw)
denv.replicates.norm <- replicates(denv.norm)
plot(denv.replicates.raw)
plot(denv.replicates.norm)

plot(denv.hyper.ana, readout.thresh=2, sig.thresh=.2)
plot(denv.priorit)

