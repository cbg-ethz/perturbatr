library(data.table, quietly=T)
library(dplyr, quietly=T)
library(dtplyr, quietly=T)
library(tidyr, quietly=T)
library(reshape2, quietly=T)
library(svd)

rm(list=setdiff(ls(), "rnai.screen.raw"))
load("inst/extdata/rnai.screen.raw.rda")

hcv.raw        <- filter(rnai.screen.raw, Virus=="HCV", Screen=="Genome")

hcv.replicates.raw <- replicates(hcv.raw)
plot(hcv.replicates.raw)

hcv.rand.lpla <- dplyr::filter(hcv.raw, Replicate==1) %>%
    dplyr::filter(Plate==3)
hcv.plates.raw <- plates(hcv.rand.lpla)
svd:::.multiplot(
  plot(hcv.plates.raw[[1]], show.gene.names=T),
  plot(hcv.plates.raw[[2]], show.gene.names=T))


hcv.quality.raw <- quality(hcv.raw)
plot(hcv.quality.raw)

hcv.norm       <- svd::preprocess(hcv.raw, normalize=c(  "b.score",  "robust-z.score"),
                                   method="mean",
                                   z.score.level="plate",
                                   rm.cytotoxic=NULL,
                                   drop=T, normalize.viability=F)

hcv.quality.norm <- quality(hcv.norm)
plot(hcv.quality.norm)

hcv.replicates.norm <- replicates(hcv.norm)
svd::cor(hcv.replicates.norm) %>% data.table::melt() %>% summarize(m=mean(value))


hcv.hyper.ana <- hyperstatistic(hcv.norm, level="gene")

hcv.t.ana     <- tstatistic(hcv.norm)
plot(hcv.t.ana)

hcv.priorit.h <- prioritize(hcv.hyper.ana, hit.ratio=.6, readout.thresh=2)
hcv.priorit.t <- prioritize(hcv.t.ana, hit.ratio=.6, readout.thresh=2)

hist(hcv.raw)
hist(hcv.norm)

hcv.plates.raw <- plates(hcv.raw)
hcv.plates.norm <- plates(hcv.norm)
sam <- sample(length(hcv.plates.raw), 1)
plot(hcv.plates.raw, count=3,  show.controls=T)
plot(hcv.plates.norm , count=3, show.controls=T)
svd:::.multiplot(
  plot(hcv.plates.raw[[sam]], show.gene.names=T),
  plot(hcv.plates.norm[[sam]], show.gene.names=T))

hcv.replicates.raw <- replicates(hcv.raw)
hcv.replicates.norm <- replicates(hcv.norm)
plot(hcv.replicates.raw)
plot(hcv.replicates.norm)

plot(hcv.hyper.ana, readout.thresh=2, sig.thresh=.2)
plot(hcv.priorit)

