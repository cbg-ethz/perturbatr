#################################
#
# read CHIKV kinase data
# this is a CONST FILE reader and should actually not be used since the data is also in data/chikv.kinase.screen.Rd

require(data.table)
require(dplyr)
require(reshape2)
require(tidyr)
require(gdata)
rm(list=ls())
options(warn=1)
assign("last.warning", NULL, envir=baseenv())


chikv.dir   <- '/Users/simondi/PHD/data/data/sysvirdrug/screening_data/chikv_kinase_screen/'
readout.dir <- paste(chikv.dir, 'CHIKV_kinase_screen_results/', sep='')
layout.file <- paste(chikv.dir, 'kinases_dharmacon_siRNA_info.xls', sep='')

read.layout <- function(lib)
{
  layout <- read.xls(lib)[c(-1, -2), -10]
  colnames(layout) <-
    c('Plate', 'Well', 'PoolID', 'siRNAID',
      'GeneSymbol', 'GeneID', 'RefSeq', 'GI', 'Sequence')
  layout$Plate <-
    as.numeric(gsub("[^0-9]+", '', layout$Plate))
  layout <-
    as.data.table(layout)
  layout <- layout  %>%
    mutate(RowIdx=match(tolower(gsub("[0-9]", "", Well)), letters),
           ColIdx=as.numeric(gsub("[^0-9]", "", Well)))
  layout
}

# parse the replicate name from the file
replicate.name <- function(f)
{
  basename <-  gsub('.xls', '',  basename(f))
  replicate <- as.integer(gsub("[^0-9]", '', basename))
  replicate
}

read.gfp <- function(f, sheet=4, cols=c(2,3,4,6,7,8))
{
  readout <- as.data.table(read.xls(f, sheet=sheet)[c(-1), cols])
  colnames(readout) <- c("Plate", "Well", "GeneSymbol", paste("Tri", 1:3, sep=""))
  readout$Plate <-  as.numeric(gsub("[^0-9]+", '', readout$Plate))
  readout <- readout  %>%
    filter(!is.na(Plate)) %>%
    mutate(RowIdx=match(tolower(gsub("[0-9]", "", Well)), letters),
           ColIdx=as.numeric(gsub("[^0-9]", "", Well)))
}

# read well gene mappings
read.target.overview<- function(f)
{
  target <- as.data.table(read.xls(f)[c(-1,-2), c(-2, -8)])
  colnames(target) <- c("Plate", "Well", "PoolID", "GeneSymbol", "LocusID", "RefSeq")
  target$Plate <-  as.numeric(gsub("[^0-9]+", '', target$Plate))
  target <- target  %>%
    mutate(RowIdx=match(tolower(gsub("[0-9]", "", Well)), letters),
           ColIdx=as.numeric(gsub("[^0-9]", "", Well)))
  target
}

get.readout.table <- function(readout.dir, sheet=4, cols=c(2,3,4,6,7,8))
{
  files <- list.files(readout.dir, full.names=T)
  readout.table   <- data.table()
  for (i in  seq(length(files)))
  {
    f <- files[i]
    replicate <- replicate.name(f)
    target <- read.target.overview(f)
    readout <- read.gfp(f, sheet=sheet, cols)
    suppressWarnings(joined <- as.data.table(
      full_join(readout, target, by=c('RowIdx', 'ColIdx', 'Plate',
                                      'Well', "GeneSymbol"))))
    joined[, 'Replicate':= as.integer(replicate)]
    readout.table <- rbindlist(list(readout.table, joined))
  }
  readout.table
}

parse.uniform <- function(readout.table)
{
  readout.table[, 'Control' := 0]
  setDT(readout.table)[Control == 0, Info := "regsample"]

  setDT(readout.table)[startsWith(GeneSymbol, "NTP"), GeneSymbol := "Scrambled"]
  setDT(readout.table)[GeneSymbol == "Scrambled", Control := -1]
  setDT(readout.table)[GeneSymbol == "Scrambled", Info:="negcontrol"]

  setDT(readout.table)[startsWith(GeneSymbol, "GAPDH"), Control := -1]
  setDT(readout.table)[startsWith(GeneSymbol, "GAPDH"), Info := "negcontrol"]

  setDT(readout.table)[startsWith(GeneSymbol, "no si"), GeneSymbol := "empty"]
  setDT(readout.table)[startsWith(GeneSymbol, "empty"), Control := as.integer(NA)]
  setDT(readout.table)[startsWith(GeneSymbol, "empty"), Info := as.character(NA)]

  setDT(readout.table)[startsWith(GeneSymbol, "buffer"), Control := -1]
  setDT(readout.table)[startsWith(GeneSymbol, "buffer"), Info := "negcontrol"]

  setDT(readout.table)[startsWith(GeneSymbol, "EAV"), Control := -1]
  setDT(readout.table)[startsWith(GeneSymbol, "EAV"), Info := "negcontrol"]
  setDT(readout.table)[startsWith(GeneSymbol, "EAV"), GeneSymbol := "EAV-nsp7"]

  readout.table
}

read.raw.data <- function(readout.data, sheet=14)
{
  files <- list.files(readout.data, full.names=T)
  full.raw.chikv.table <- data.table()
  for (i in  seq(length(files)))
  {
    xls <- files[i]
    raw.chikv.data <- as.data.table(read.xls(xls,sheet=sheet))
    raw.chikv.data[,
                  which(
                    unlist(
                      lapply(raw.chikv.data,
                             function(x) all(is.na(x) | x == 0))))
                  :=NULL,
                  with=F]
    raw.chikv.data <- raw.chikv.data[,rowsum := rowSums(.SD), .SDcols=-1] %>%
      dplyr::filter(rowsum != 0) %>%
      dplyr::select(-rowsum)
    colnames(raw.chikv.data) <- c("Plate",
                                 rep("Triplicate1", 12),
                                 rep("Triplicate2", 12),
                                 rep("Triplicate3", 12))
    plates <- as.integer(gsub("[^0-9]" ,"", raw.chikv.data$Plate))
    ind <- which(!is.na(plates))
    plates <- rep(plates[ind], times=diff(c(ind, length(plates)+1)))
    raw.chikv.data[,Plate := plates]
    formatted.chikv.data <- raw.chikv.data[,.SD,.SDcols=1:13]
    a <- cbind(1, raw.chikv.data[,.SD,.SDcols=1:13])
    b <- cbind(2, raw.chikv.data[,.SD,.SDcols=c(1,14:25)])
    c <- cbind(3, raw.chikv.data[,.SD,.SDcols=c(1,26:37)])
    colnames(a) <- colnames(b) <- colnames(c) <-  c("Triplicate", "Plate", 1:12)
    formatted.chikv.data <-  rbind(a,b,c)
    formatted.chikv.data <- formatted.chikv.data %>%
      group_by(Triplicate, Plate) %>%
      .[,RowIdx:=1:8]
    formatted.chikv.data <- as.data.frame(dplyr::ungroup(formatted.chikv.data))
    final.data <- data.table()
    n.row <- nrow(formatted.chikv.data)
    n.col <- ncol(formatted.chikv.data)
    for (ridx in seq(n.row))
    {
      rep <- formatted.chikv.data$Triplicate[ridx]
      plate <- formatted.chikv.data$Plate[ridx]
      rowidx <- formatted.chikv.data$RowIdx[ridx]
      for (col in 3:14)
      {
        if (ridx == 1 & col == 3){
          final.data <- as.data.table(list(Replicate=i, Triplicate=rep,
                                           Plate=plate, RowIdx=rowidx,
                                           ColIdx=col-2,
                                           Readout=formatted.chikv.data[ridx,col]))
        }
        else {
          final.data <- rbindlist(list(final.data,
                                       list(Replicate=i, rep, plate,
                                            rowidx, col-2,
                                            formatted.chikv.data[ridx,col])))
        }
      }
    }
    final.data <- final.data %>%
      group_by(Plate, RowIdx, ColIdx) %>%
      spread(Triplicate, Readout)
    final.data <- ungroup(final.data)
    colnames(final.data) <-
      c("Replicate", "Plate", "RowIdx", "ColIdx", paste("Tri", 1:3, sep=""))
    if (i == 1) {
      full.raw.chikv.table <- final.data
    } else {
      full.raw.chikv.table <- rbindlist(list(full.raw.chikv.table, final.data))
    }
  }

  full.raw.chikv.table
}

add.raw.data.to.readout <- function(readout.table, raw.data.table)
{
  n.row <- nrow(raw.data.table)
  n.row.prev <- nrow(readout.table)
  sucs <- fls <- 0
  amb.wells <- c()
  for (i in seq(nrow(raw.data.table)))
  {

    repl <- raw.data.table$Replicate[i]
    plate <- raw.data.table$Plate[i]
    rowidx <- raw.data.table$RowIdx[i]
    colidx <- raw.data.table$ColIdx[i]
    raw.readout <- c(raw.data.table$Tri1[i], raw.data.table$Tri2[i],  raw.data.table$Tri3[i])
    val <- filter(readout.table, Plate==plate, Replicate==repl, RowIdx==rowidx, ColIdx==colidx)
    if (!is.na(val) && nrow(val) != 0 )
    {
      if (any(c(val$Tri1, val$Tri2, val$Tri3) != raw.readout))
      {
        cat("Wrong values", i, "\n" )
        stop("There is an error with data integration!")
      } else {
        sucs <- sucs + 1
      }
    } else {
      well <- paste(toupper(letters[rowidx]), sprintf("%02d", colidx), sep="")
      amb.wells <- c(amb.wells, well)
      fls <- fls + 1
      if (all(c("PoolID", "LocusID") %in% colnames(readout.table))) {
        readout.table <- rbindlist(
          list(readout.table, list(plate, well, "undef",
                                   raw.readout[1], raw.readout[2], raw.readout[3],
                                   rowidx, colidx,
                                   as.character(NA), as.character(NA), as.character(NA),
                                   repl, as.character(NA), as.character(NA) )))
      } else {
        readout.table <- rbindlist(
          list(readout.table, list(repl, plate, well, "undef",
                                   rowidx, colidx,
                                   raw.readout[1], raw.readout[2], raw.readout[3])))
      }

    }
  }
  cat("Success", sucs, "\n")
  cat("Fails", fls, "\n")
  cat(paste(unique(amb.wells), collapse=", "), "\n")
  if (sucs != n.row.prev) stop("You made a parsing error, brah!")
  if (fls != n.row - n.row.prev) stop("You made a parsing error, fool!")
  cat("Sweeeeet!\n")
  readout.table
}

read.primary.readout <- function(readout.dir)
{

  readout.table  <- get.readout.table(readout.dir)
  readout.table <- parse.uniform(readout.table)
  raw.readout.table <- read.raw.data(readout.dir)
  readout.table <- add.raw.data.to.readout(readout.table, raw.readout.table)
  readout.table

}

join.readout.layout <-
function(readout.table, layout.table)
{
  merged.layout.table <-
    dplyr::select(layout.table, Plate, Well, PoolID, siRNAID,
           GeneSymbol, GeneID, RefSeq,
           RowIdx, ColIdx, Sequence) %>%
    dplyr::group_by( Plate, Well, GeneSymbol, PoolID, RefSeq, GeneID, RowIdx, ColIdx) %>%
    dplyr::summarize(Sequences=paste(Sequence, collapse=","), siRNAIDs=paste(siRNAID, collapse=",")) %>%
    ungroup %>%
    dplyr::mutate(RowIdx=match(tolower(gsub("[0-9]", "", Well)), letters),
           ColIdx=as.numeric(gsub("[^0-9]", "", Well)) - 1) %>% # subtract 1 due to rearrangement after transfecter
    dplyr::select(GeneSymbol, PoolID, RefSeq, Sequences, siRNAIDs, RowIdx, ColIdx)
  final.table <- suppressWarnings(
    as.data.table(full_join(readout.table,
                            merged.layout.table,
                            by=c("RowIdx", "ColIdx", "PoolID"))))
  setDT(final.table)[GeneSymbol.x == GeneSymbol.y,  GeneSymbol.y:=as.character(NA)]
  setDT(final.table)[as.character(RefSeq.x) == as.character(RefSeq.y),  RefSeq.y:=as.character(NA)]
  final.table <- filter(final.table, !is.na(Plate) & !is.na(Well))
  colnames(final.table)[which(startsWith(colnames(final.table), "Tri"))] <-
    paste("Readout", 1:3, sep="")
  colnames(final.table)[which(startsWith(colnames(final.table), "Gene"))] <-
    c("GeneSymbol","AltGeneSymbol")
  colnames(final.table)[which(startsWith(colnames(final.table), "RefSeq"))] <-
    c("RefSeq","AltRefSeq")
  final.table
}

read.viability <- function(readout.dir)
{
  viability.table  <- get.readout.table(readout.dir, sheet=3, cols=c(1,2,3,5,6,7))
  viability.table  <- filter(viability.table, !is.na(Plate) & !is.na(Well)) %>%
    dplyr::select(Replicate, Plate, Well, GeneSymbol,
           RowIdx, ColIdx, Tri1, Tri2, Tri3)
  viability.table <- parse.uniform(viability.table)
  viability.table <- dplyr::select(viability.table, Replicate, Plate, Well,
                            GeneSymbol, RowIdx, ColIdx,
                            Tri1, Tri2, Tri3)
  raw.readout.table <- read.raw.data(readout.dir, sheet=13)
  viability.table <- add.raw.data.to.readout(viability.table, raw.readout.table)
  colnames(viability.table)[which(startsWith(colnames(viability.table), "Tri"))] <-
    paste("Viability", 1:3, sep="")
  viability.table
}

join.readout.viability <-
function(readout.table, viability.table)
{
    final.table <- as.data.table(full_join(joined.primary.layout.table,
                                           viability.table,
                                           by=c("Replicate", "Plate", "Well", "RowIdx",
                                                "ColIdx", "GeneSymbol")))
    final.table[, 'Library':= 'DharmaconSMARTPool']
    setDT(final.table)[startsWith(GeneSymbol, "buffer"), GeneSymbol:="buffer"]
    setcolorder(final.table,
                c("Replicate", "Well", "RowIdx", "ColIdx", "Plate", "GeneSymbol",
                  paste("Readout",1:3,sep=""), paste("Viability",1:3,sep=""),
                  "Control", "Info", "Library", "RefSeq", "AltGeneSymbol", "AltRefSeq",
                  "LocusID",  "PoolID",
                  "Sequences", "siRNAIDs"))
    final.table <- final.table[order(Well, Plate, GeneSymbol)]
    final.table
}

post.process <- function(screen)
{
  setDT(screen)[, Control := as.numeric(Control)]

  setDT(screen)[RowIdx == 1 & ColIdx==1, GeneSymbol := "undef"]
  setDT(screen)[RowIdx == 2 & ColIdx==1, GeneSymbol := "buffer"]
  setDT(screen)[RowIdx == 3 & ColIdx==1, GeneSymbol := "Scrambled"]
  setDT(screen)[RowIdx == 4 & ColIdx==1, GeneSymbol := "GAPDH"]
  setDT(screen)[RowIdx == 5 & ColIdx==1, GeneSymbol := "buffer"]
  setDT(screen)[RowIdx == 6 & ColIdx==1, GeneSymbol := "Scrambled"]
  setDT(screen)[RowIdx == 7 & ColIdx==1, GeneSymbol := "GAPDH"]
  setDT(screen)[RowIdx == 8 & ColIdx==1, GeneSymbol := "buffer"]

  setDT(screen)[startsWith(GeneSymbol, "NTP"), GeneSymbol := "Scrambled"]
  setDT(screen)[GeneSymbol == "Scrambled", Control := -1]
  setDT(screen)[GeneSymbol == "Scrambled", Info:="negcontrol"]

  setDT(screen)[startsWith(GeneSymbol, "GAPDH"), Control := -1]
  setDT(screen)[startsWith(GeneSymbol, "GAPDH"), Info := "negcontrol"]

  setDT(screen)[startsWith(GeneSymbol, "no si"), GeneSymbol := "empty"]
  setDT(screen)[startsWith(GeneSymbol, "empty"), Control := as.integer(NA)]
  setDT(screen)[startsWith(GeneSymbol, "empty"), Info := as.character(NA)]

  setDT(screen)[startsWith(GeneSymbol, "buffer"), Control := -1]
  setDT(screen)[startsWith(GeneSymbol, "buffer"), Info := "negcontrol"]

  setDT(screen)[startsWith(GeneSymbol, "EAV"), Control := -1]
  setDT(screen)[startsWith(GeneSymbol, "EAV"), Info := "negcontrol"]
  setDT(screen)[startsWith(GeneSymbol, "EAV"), GeneSymbol := "EAV-nsp7"]

  setDT(screen)[startsWith(GeneSymbol, "undef"), GeneSymbol := NA_character_]
  setDT(screen)[startsWith(GeneSymbol, "empty"), GeneSymbol := NA_character_]

  screen
}

layout.table    <- read.layout(layout.file)
primary.readout.table   <- read.primary.readout(readout.dir)
joined.primary.layout.table  <-
  join.readout.layout(primary.readout.table,
                      layout.table)
viability.table <- read.viability(readout.dir)
final.table <- join.readout.viability(
                joined.primary.layout.table,
                viability.table)

final.table <- post.process(final.table)

chikv.kinase.screen <- final.table
save(file="inst/extdata/chikv.kinase.screen.rda", chikv.kinase.screen, compress="bzip2")
write.table(chikv.kinase.screen,
           file="~/PHD/data/data/sysvirdrug/integrated_data_files/chikv_kinase_screen.tsv",
           sep='\t', quote=F, row.names=F)

#delete all and load data
rm(list=ls())
load('inst/extdata/chikv.kinase.screen.rda')
