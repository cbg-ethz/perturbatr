##################
# Final data integration!
##################

library(data.table)
library(dplyr)
library(dtplyr)
library(reshape2)
library(gdata)
library(tidyr)
rm(list=ls())
options(warn=1)
assign("last.warning", NULL, envir=baseenv())

# cvb not complete yet
#load("data/cvb.replicon.screen.Rd")

join.denv.and.hcv.kinase.screen <- function()
{
  load("~/PHD/data/data/sysvirdrug/integrated_data_files/hcv.kinase.screen.rda")
  load("~/PHD/data/data/sysvirdrug/integrated_data_files/denv.kinase.screen.rda")
  denv.kinase.screen[, 'Virus' := 'DENV']
  denv.kinase.screen[, 'Screen' := 'Kinome']
  hcv.kinase.screen[, 'Virus' := 'HCV']
  hcv.kinase.screen[, 'Screen' := 'Kinome']

  kinase.screen <- rbindlist(list(hcv.kinase.screen, denv.kinase.screen))
  first.cols <- c("Virus", "Replicate", "Plate", "WellNo", "RowIdx", "ColIdx")
  res.cols <- setdiff(colnames(kinase.screen), first.cols)
  setcolorder(kinase.screen, c(first.cols, res.cols))
  rm(hcv.kinase.screen)
  rm(denv.kinase.screen)

  kinase.screen
}

join.denv.and.hcv.druggable.screen <- function()
{
  load("~/PHD/data/data/sysvirdrug/integrated_data_files/denv.druggable.genome.screen.rda")
  load("~/PHD/data/data/sysvirdrug/integrated_data_files/hcv.druggable.genome.screen.rda")
  denv.druggable.screen[, 'Virus' := 'DENV']
  denv.druggable.screen[, 'Screen' := 'Genome']
  hcv.druggable.screen[, 'Virus' := 'HCV']
  hcv.druggable.screen[, 'Screen' := 'Genome']

  druggable.genome.screen <- rbindlist(list(hcv.druggable.screen, denv.druggable.screen))
  first.cols <- c("Virus", "Replicate", "Plate", "Well", "RowIdx", "ColIdx")
  res.cols <- setdiff(colnames(druggable.genome.screen), first.cols)
  setcolorder(druggable.genome.screen, c(first.cols, res.cols))
  rm(hcv.druggable.screen)
  rm(denv.druggable.screen)
  setDT(druggable.genome.screen)[,InfectionType := gsub("k", "c", InfectionType)]
  setDT(druggable.genome.screen)[,Library := "Ambion"]

  druggable.genome.screen
}

join.denv.and.hcv.screens <- function()
{
  # TODO are well indexes correct?
  kinase.screen <- join.denv.and.hcv.kinase.screen()
  setDT(kinase.screen)[,ReadoutType := "GFP"]
  setDT(kinase.screen)[,InfectionType := "Infection/Reinfection"]
  kinase.screen <- dplyr::select(kinase.screen, -Sequence, -WellNo, -RefSeq)
  genome.screen <- join.denv.and.hcv.druggable.screen()
  colnames(kinase.screen)[startsWith(colnames(kinase.screen), "MeanCyt" )] <- "Readout"
  colnames(kinase.screen)[startsWith(colnames(kinase.screen), "siRNA" )] <- "siRNAID"
  genome.screen <- dplyr::select(genome.screen, -GeneID, -PlateName, -Sequence, -Well, -RefSeq)
  setDT(genome.screen)[, NumCells:=as.integer(NA) ]
  setcolorder(genome.screen, colnames(kinase.screen))
  screen <- rbindlist(list(kinase.screen, genome.screen))
  screen[, 'Cell' := 'HuH7.5']
  screen
}

rearrange.triplicates <- function(obj, ...)
{
  max.rep <- max(obj$Replicate)
  max.tri <- 3
  new.rep <- function(rep, read)
  {
    read <-as.integer(sub("[[:alpha:]]+", "", read))
    ((rep - 1) * max.tri) + read
  }
  s <- tidyr::gather(obj, ReadoutClass, Readout, c(Readout1, Readout2, Readout3,
                                    Viability1, Viability2, Viability3)) %>%
    as.data.table
  s <- dplyr::mutate(s, Replicate=new.rep(Replicate, ReadoutClass)) %>% as.data.table
  s <- dplyr::mutate(s, ReadoutClass=sub("[[:digit:]]", "", ReadoutClass)) %>% as.data.table %>%
    .[order(Replicate, Plate, RowIdx, ColIdx)]
  s
}

join.sars.and.chikv.screens <- function()
{
  load("~/PHD/data/data/sysvirdrug/integrated_data_files/chikv.kinase.screen.rda")
  # TODO sars kinase not done complete yet
  load("~/PHD/data/data/sysvirdrug/integrated_data_files/sars.kinase.screen.rda")

  chikv.kinase.screen[, 'Virus' := 'CHIKV']
  chikv.kinase.screen[, 'Screen' := 'Kinome']
  chikv.kinase.screen[, 'Cell' := 'MRC5']
  sars.kinase.screen[, 'Virus' := 'SARS']
  sars.kinase.screen[, 'Screen' := 'Kinome']
  sars.kinase.screen[, 'Cell' := 'ACE2']

  chikv.kinase.screen <- dplyr::select(chikv.kinase.screen,
                                       -RefSeq, -Well, -AltGeneSymbol,
                                       -AltRefSeq, -LocusID, -PoolID, -Sequences)
  chikv.kinase.screen <- rearrange.triplicates(chikv.kinase.screen)
  sars.kinase.screen  <- dplyr::select(sars.kinase.screen,
                                       -RefSeq, -Well, -AltGeneSymbol,
                                       -AltRefSeq, -LocusID, -PoolID, -Sequences)
  sars.kinase.screen <- rearrange.triplicates(sars.kinase.screen)

  setcolorder(sars.kinase.screen, colnames(chikv.kinase.screen))
  kinase.screen <- rbindlist(list(sars.kinase.screen, chikv.kinase.screen))
  rm(chikv.kinase.screen)
  rm(sars.kinase.screen)
  setDT(kinase.screen)[,ReadoutType := "GFP"]
  setDT(kinase.screen)[,InfectionType := "Infection/Reinfection"]
  imp.cols <- c("Virus", "Replicate", "Plate", "RowIdx", "ColIdx")
  res.cols <- setdiff(colnames(kinase.screen), imp.cols)
  setcolorder(kinase.screen, c(imp.cols, res.cols))
  kinase.screen
}

join.sc.and.hd.screens <- function()
{
  sc.screen <- join.sars.and.chikv.screens()
  hd.screen <- join.denv.and.hcv.screens()

  setDT(sc.screen)[, "NumCells" := as.integer(NA)]
  setDT(hd.screen)[, "ReadoutClass" := "Readout"]
  setDT(hd.screen)[ReadoutType == "F-Luciferase", "ReadoutClass" := "Viability"]
  setDT(hd.screen)[ReadoutType == "F-Luciferase", ReadoutType := "Luciferase"]
  setDT(hd.screen)[ReadoutType == "R-Luciferase", ReadoutType := "Luciferase"]
  colnames(hd.screen)[which(startsWith(colnames(hd.screen), "siR"))] <- "siRNAIDs"
  col.nms.sc <- colnames(sc.screen)
  col.nms.hd <- colnames(hd.screen)
  setcolorder(hd.screen, col.nms.sc)
  screen <- rbindlist(list(hd.screen, sc.screen))
  imp.cols <- c("Virus", "Replicate", "Plate",
                "RowIdx", "ColIdx", "GeneSymbol",
                colnames(screen)[startsWith(colnames(screen), "Readout")],
                colnames(screen)[startsWith(colnames(screen), "Viability")])
  res.cols <- setdiff(colnames(screen), imp.cols)
  setcolorder(screen, c(imp.cols, res.cols))
  screen
}

integrate <- function()
{
  screen <- join.sc.and.hd.screens()
  imp.cols <- c("Virus", "Replicate", "Plate",
                "RowIdx", "ColIdx", "GeneSymbol",
                "ReadoutClass", "Readout")

  setDT(screen)[,GeneSymbol := as.character(GeneSymbol)]
  setDT(screen)[,ReadoutClass := as.character(ReadoutClass)]
  setDT(screen)[,siRNAIDs := as.character(siRNAIDs)]
  setDT(screen)[,Plate := as.integer(Plate)]
  setDT(screen)[,Replicate := as.integer(Replicate)]
  setDT(screen)[,NumCells := as.integer(NumCells)]
  setDT(screen)[,Control := as.integer(Control)]
  setDT(screen)[,RowIdx := as.integer(RowIdx)]
  setDT(screen)[,ColIdx := as.integer(ColIdx)]
  res.cols <- setdiff(colnames(screen), imp.cols)
  setcolorder(screen, c(imp.cols, res.cols))
  setDT(screen)[InfectionType=="Infection", InfectionType := "I"]
  setDT(screen)[InfectionType=="Reinfection", InfectionType := "R"]
  setDT(screen)[InfectionType=="Infection/Reinfection", InfectionType := "IR"]
  setDT(screen)[,Meta := paste(Virus, "Replicate", Replicate, "Plate", Plate, ReadoutClass, ReadoutType, InfectionType, Screen, sep="_")]

  setDT(screen)[GeneSymbol == "Trim5alpha", GeneSymbol := "TRIM5"]
  setDT(screen)[GeneSymbol == "p53", GeneSymbol := "TP53"]
  setDT(screen)[GeneSymbol == "LEDGF-p75", GeneSymbol := "PSIP1"]
  setDT(screen)[GeneSymbol == "PP1A", GeneSymbol := "PPP1CA"]
  setDT(screen)[GeneSymbol == "PP1B", GeneSymbol := "PPP1CB"]
  setDT(screen)[GeneSymbol == "PP1C", GeneSymbol := "PPP1CC"]
  setDT(screen)[GeneSymbol == "CKI-delta", GeneSymbol := "CSNK1D"]
  setDT(screen)[GeneSymbol == "CKI-alpha", GeneSymbol := "CSNK1A1"]
  setDT(screen)[GeneSymbol == "CKI-alpha", GeneSymbol := "CSNK1A1"]
  setDT(screen)[GeneSymbol == "bCop", GeneSymbol := "COPB"]

  setDT(screen)[                 , Design := "pooled"]
  setDT(screen)[Library=="Ambion", Design := "single"]

  setDT(screen)[Info=="<NA>"    , Info := NA_character_]
  setDT(screen)[siRNAIDs=="<NA>", siRNAIDs := NA_character_]
  setDT(screen)[GeneSymbol=="<NA>", GeneSymbol := NA_character_]
  invisible(as.data.table(screen))
}

map.entrez <- function(screen)
{
  load.hugo.mappings <-
  function ()
  {
    hugo2entrez <-
      as.data.table(read.table("/Users/simondi/PHD/data/data/sysvirdrug/mappings/hugo2entrez_r_mappings_integration.tsv",
                               sep="\t", fill=T, header=F))
    colnames(hugo2entrez) <- c("GeneSymbol", "Entrez", "Source")
    setDT(hugo2entrez)[ ,Source := NULL]
    hugo2entrez$GeneSymbol <- as.character(hugo2entrez$GeneSymbol)
    hugo2entrez
  }

  #' Join the hugo2Entrez mapping with the RNAi data.table
  #'
  #' @noRd
  #' @importFrom dplyr left_join
  #' @importFrom data.table setDT
  map <- function(dat)
  {
    hugo2entrez <- load.hugo.mappings()
    cols <- colnames(dat)
    dat <- dplyr::left_join(dat, hugo2entrez, by="GeneSymbol")
    setcolorder(dat, c(cols, "Entrez"))
    setDT(dat)[GeneSymbol == "<NA>", GeneSymbol := NA_character_]
    setorder(dat, Virus, Replicate, Plate, RowIdx, ColIdx)
    dat
  }
  screen <- map(screen)
  screen
}

finalize <- function(screen)
{
  setDT(screen)[GeneSymbol == "Scrambled", Control := -1]
  setDT(screen)[GeneSymbol == "Scrambled", Info:="negcontrol"]

  screen
}

write.screen <- function(screen)
{
  rnai.screen.raw <- screen
  save(file="./inst/extdata/rnai.screen.raw.rda",
       rnai.screen.raw,
       compress="bzip2")
  write.table(rnai.screen.raw,
              file="~/PHD/data/data/sysvirdrug/integrated_data_files/rnai.screen.raw.tsv",
              sep='\t',quote=F,row.names=F)
}

screen <- integrate()
screen <- map.entrez(screen)
screen <- finalize(screen)
screen <- as.data.table(screen)  %>% dplyr::select(-Meta, -Info)
class(screen) <- c("svd.raw", class(screen))


#delete all and load data
write.screen(screen)
rm(list=ls())
load('./inst/extdata/rnai.screen.raw.rda')
