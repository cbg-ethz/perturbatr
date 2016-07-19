#################################
#
# read DENV druggable genome data
# this is a CONST FILE reader and should actually not be used since the data is also in data/denv.druggable.genome..screen.Rd

require(data.table)
require(dplyr)
require(reshape2)
require(tidyr)
require(gdata)
library(stringr)
rm(list=ls())
options(warn=2)
assign("last.warning", NULL, envir=baseenv())

denv.dir <- '/Users/simondi/PHD/data/data/sysvirdrug/screening_data/dengue_druggable_genome_screen/'
lib.data <- paste(denv.dir,  'ambion_library_layout/', sep='')
denv.data.dir <- paste(denv.dir,  'DEN_druggable_genome_results/', sep='')

parse.uniform.layout <- function(full.layout.table)
{

  colnames(full.layout.table) <-
    c("Plate", "Well", "GeneSymbol", "GeneID", "siRNAID", "RefSeq",
      "Sequence", "IsValidated", "negControl", "posControl")

  full.layout.table <- filter(full.layout.table, GeneSymbol != "")

  setDT(full.layout.table)[, GeneSymbol := as.character(GeneSymbol)]
  setDT(full.layout.table)[, RefSeq := as.character(RefSeq)]
  setDT(full.layout.table)[RefSeq=="n/a", RefSeq := as.character(NA)]
  setDT(full.layout.table)[, Sequence := as.character(Sequence)]
  setDT(full.layout.table)[, Well := as.character(trim(Well))]
  setDT(full.layout.table)[startsWith(GeneSymbol, "Empty") , GeneSymbol := as.character("empty")]
  setDT(full.layout.table)[RefSeq == "" , RefSeq := as.character(NA)]
  setDT(full.layout.table)[ , Library := as.character("Ambion")]
  setDT(full.layout.table)[ , Control := 0]
  setDT(full.layout.table)[ , Info := "regsample"]
  setDT(full.layout.table)[posControl == 1, Control :=  1]
  setDT(full.layout.table)[posControl == 1, Info := "poscontrol"]
  setDT(full.layout.table)[negControl == 1, Control := -1]
  setDT(full.layout.table)[negControl == 1, Info := "negcontrol"]


  full.layout.table <- full.layout.table  %>%
    mutate(RowIdx=match(tolower(gsub("[0-9]", "", Well)), letters),
           ColIdx=as.integer(gsub("[^0-9]", "", Well)))
  full.layout.table <- dplyr::select(full.layout.table, Plate, Well, RowIdx, ColIdx, GeneSymbol, GeneID,
                              siRNAID, RefSeq, Sequence, Control, Library, Info) %>%
                       dplyr::mutate(Plate=as.integer(gsub("[^0-9]", "", Plate)))


  control.genes <- filter(full.layout.table, Control != 0) %>%
    dplyr::select(GeneSymbol) %>% unique %>% unlist %>% as.vector
  pos.control.genes <- control.genes[startsWith(control.genes, "DV")]
  neg.control.genes <- setdiff(control.genes, pos.control.genes)

  setDT(full.layout.table)[GeneSymbol %in% pos.control.genes, Control:= 1]
  setDT(full.layout.table)[GeneSymbol %in% pos.control.genes, Info:="poscontrol"]
  setDT(full.layout.table)[GeneSymbol %in% neg.control.genes, Control:=-1]
  setDT(full.layout.table)[GeneSymbol %in% neg.control.genes, Info:="negcontrol"]
  full.layout.table
}

read.layout <- function(lib)
{
  files <- list.files(lib, full.names=T)
  full.layout.table <- data.table()

   for (i in seq(length(files)))
  {
    file <- files[i]
    layout.table <- as.data.table(read.xls(file))
    full.layout.table <- rbindlist(list(full.layout.table, layout.table),fill=T)
  }
  full.layout.table <- parse.uniform.layout(full.layout.table)
  setkey(full.layout.table, NULL)
  unique(full.layout.table)
}

read.single.file <- function(f)
{
  readout.table <- as.data.table(read.table(f, sep="\t", header=F, fill=T)[-1,])
  colnames(readout.table) <- c("Well", "Readout")
  readout.table <- mutate(readout.table,
                          RowIdx=match(tolower(gsub("[0-9]", "", trim(Well))), letters),
                          ColIdx=as.integer(gsub("[^0-9]", "", trim(Well))))
  readout.table
}

plate.info <- function(f)
{
    f <- gsub(".*/", "", f)
    regex <- "([:alpha:]+([:digit:]+)\\-([:digit:]+)\\_([:alpha:]+)\\_([:alpha:]{1})).[:alpha:]{3}"
    match <- str_match(f, regex)
    if (match[5] == "Reinf"){ s <- "Reinfection"}
    else if (match[5] == "Inf"){ s <- "Infection"}
    else {stop("Parsing error")}
    if (match[6] == "R"){ inf <- "R-Luciferase"}
    else if (match[6] == "F"){ inf <- "F-Luciferase"}
    else {stop("Parsing error")}
    list(plate.name=match[2],
         replicate=as.integer(match[3]),
         plate.idx=as.integer(match[4]),
         infection.type=s,
         readout.type=inf)
}

read.primary.readout <- function(denv.data.dir)
{
  files <- list.files(denv.data.dir, recursive=T, full.names=T)
  full.readout.table <- data.table()
  for (i in seq(length(files)))
  {
    file <- files[i]
    plate.info <- plate.info(file)
    readout.table <- read.single.file(file)
    setDT(readout.table)[,PlateName := as.character(plate.info$plate.name)]
    setDT(readout.table)[,Plate := plate.info$plate.idx]
    setDT(readout.table)[,Well := trim(Well)]
    setDT(readout.table)[,ReadoutType := as.character(plate.info$readout.type)]
    setDT(readout.table)[,InfectionType := as.character(plate.info$infection.type)]
    setDT(readout.table)[,Replicate := plate.info$replicate]
    full.readout.table <- rbindlist(list(full.readout.table , readout.table))
  }
  setDT(full.readout.table)[, Readout := as.numeric(as.character(Readout))]
  setDT(full.readout.table)[, Well := as.character(Well)]

  full.readout.table
}

join.tables <- function(readout.table, layout.table)
{

  gathered.readout.table <- as.data.table(
    full_join(readout.table, layout.table, by=c("RowIdx", "ColIdx", "Plate", "Well")))
  if (any(dim(readout.table)[1] != dim(gathered.readout.table)[1]))
  {
    stop("Join did not go well! Wrong dimensions")
  }
  gathered.readout.table
}

post.process <- function(obj, ...)
{
  setDT(obj)[startsWith(GeneSymbol, "empty"), GeneSymbol := NA_character_]
  setDT(obj)[startsWith(siRNAID, "empty"), siRNAID := NA_character_]

  setDT(obj)[startsWith(GeneSymbol, "DV-"), Control := 1]
  setDT(obj)[startsWith(GeneSymbol, "DV-"), Info := "poscontrol"]

  setDT(obj)[startsWith(GeneSymbol, "HCV"), Control := -1]
  setDT(obj)[startsWith(GeneSymbol, "HCV"), Info := "negcontrol"]

  setDT(obj)[startsWith(GeneSymbol, "Scrambled"), Control := -1]
  setDT(obj)[startsWith(GeneSymbol, "Scrambled"), Info := "negcontrol"]

  setDT(obj)[startsWith(GeneSymbol, "GFP"), Control := -1]
  setDT(obj)[startsWith(GeneSymbol, "GFP"), Info := "negcontrol"]

  setDT(obj)[startsWith(GeneSymbol, "CD4"), Control := 0]
  setDT(obj)[startsWith(GeneSymbol, "CD4"), Info := "regsample"]

  setDT(obj)[startsWith(GeneSymbol, "CD81"), Control := 0]
  setDT(obj)[startsWith(GeneSymbol, "CD81"), Info := "regsample"]

  setDT(obj)[startsWith(GeneSymbol, "GFPT1"), Control := 0]
  setDT(obj)[startsWith(GeneSymbol, "GFPT1"), Info := "regsample"]

  setDT(obj)[startsWith(GeneSymbol, "GFPT2"), Control := 0]
  setDT(obj)[startsWith(GeneSymbol, "GFPT2"), Info := "regsample"]

  setDT(obj)[startsWith(GeneSymbol, "TSG101"), Control := 0]
  setDT(obj)[startsWith(GeneSymbol, "TSG101"), Info := "regsample"]
}

layout.table <- read.layout(lib.data)
readout.table <- read.primary.readout(denv.data.dir)
joined.readout.table <- join.tables(readout.table, layout.table)
joined.readout.table <- post.process(joined.readout.table)


denv.druggable.screen <- joined.readout.table
save(file="inst/extdata/denv.druggable.genome.screen.rda", denv.druggable.screen, compress="bzip2")
write.table(denv.druggable.screen,
            file="~/PHD/data/data/sysvirdrug/integrated_data_files/denv_druggable_genome_screen.tsv",
            sep='\t',
            quote=F,
            row.names=F)
#delete all and load data
rm(list=ls())
load('inst/extdata/denv.druggable.genome.screen.rda')

