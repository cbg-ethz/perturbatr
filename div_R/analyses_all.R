library(data.table, quietly=T)
library(dplyr, quietly=T)
library(dtplyr, quietly=T)
library(tidyr, quietly=T)
library(reshape2, quietly=T)
library(svd)
library(RNAither)
library(microbenchmark)
load("inst/extdata/rnai.screen.raw.rda")
rm(list=setdiff(ls(), "rnai.screen.raw"))

dat.raw        <- filter(rnai.screen.raw, Virus=="HCV", Screen=="Kinome")

df <- filter(dat.raw, ReadoutClass=="Readout") %>% as.data.table %>%
  dplyr::group_by(Virus, Screen, Library,
         ReadoutClass, ReadoutType, InfectionType,
         Design, Cell,
         Replicate, Plate) %>%
  mutate(Spotnumber=((RowIdx-1)*max(ColIdx))+ColIdx)  %>% ungroup

df$Internal_GeneID <- df$GeneSymbol
df$GeneName <- df$GeneSymbol
df$SpotType <- -1
df$SigIntensity <- df$Readout
df$LabtekNb <- df$Plate
df$RowNb <- df$RowIdx
df$ColNb <- df$ColIdx
df$ScreenNb <- df$Replicate
df$NbCells <- df$NumCells
setDT(df)[Control == -1, SpotType := 0]
setDT(df)[Control == 1, SpotType := 1]
setDT(df)[Control == 0, SpotType := 2]
setDT(df)[, Control := NULL]
setDT(df)[, Replicate := NULL]
setDT(df)[, Plate := NULL]
setDT(df)[, RowIdx := NULL]
setDT(df)[, ColIdx := NULL]
setDT(df)[, NumCells := NULL]
setDT(df)[, siRNAIDs := NULL]
setDT(df)[, Entrez := NULL]
setDT(df)[, GeneSymbol := NULL]
setDT(df)[, ReadoutClass := NULL]
setDT(df)[, Virus := NULL]
setDT(df)[, ReadoutType := NULL]
setDT(df)[, InfectionType := NULL]
setDT(df)[, Library := NULL]
setDT(df)[, Screen := NULL]
setDT(df)[, Cell := NULL]
setDT(df)[, Design := NULL]
setDT(df)[, Readout := NULL]


rnaither(data=as.data.frame(df),expname="Example",excludeCellcounts="both",
         normalization=c("lowess", "bscore", "zscore"), scorethresh=1.0,
         outdir="~/Desktop/results/" )
