#################################
#
# read SARS-C kinase data
# this is a CONST FILE reader and should actually not be used since the data is also in data/sars.kinase.screen.Rd

require(data.table)
require(dplyr)
require(reshape2)
require(tidyr)
require(gdata)
rm(list=ls())
options(warn=1)
assign("last.warning", NULL, envir=baseenv())

sars.dir       <- '/Users/simondi/PHD/data/data/sysvirdrug/screening_data/sarscov_kinase_screen/'
lib.data       <- paste(sars.dir,  'kinases_dharmacon_siRNA_info.xls', sep='')
readout.data   <- paste(sars.dir,  'SARSCOV_kinase_screen_results/', sep='')
viability.data <- paste(sars.dir,  'viability_screen_results/', sep='')


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

# read meta data
read.sars.layout <- function(lib.data)
{
  read.layout(lib.data)
}

# parse the replicate name from the file
replicate.name <- function(f)
{
  basename <-  gsub('.xls', '',  basename(f))
  replicate <- as.integer(gsub("[^0-9]", '', basename))
  replicate
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

# read the gfp values
read.gfp <- function(f)
{
  readout <- as.data.table(read.xls(f, sheet=2)[c(-1),c(1,2,3,5,6,7)])
  colnames(readout) <- c("Plate", "Well", "GeneSymbol", paste("Tri", 1:3, sep=""))
  readout$Plate <-  as.numeric(gsub("[^0-9]+", '', readout$Plate))
  readout <- readout  %>%
    mutate(RowIdx=match(tolower(gsub("[0-9]", "", Well)), letters),
           ColIdx=as.numeric(gsub("[^0-9]", "", Well)))
  readout
}

get.readout.table <- function(readout.data)
{
  files <- list.files(readout.data, full.names=T)
  readout.table <- data.table()
  for (i in  seq(length(files)))
  {
    f <- files[i]
    replicate <- replicate.name(f)
    target <- read.target.overview(f)
    readout <- read.gfp(f)
    # join the two tables to a biiiiiig table
    suppressWarnings(joined <- as.data.table(
      full_join(readout, target, by=c('RowIdx', 'ColIdx', 'Plate', 'Well', "GeneSymbol"))))
    joined[, 'Replicate':= as.integer(replicate)]
    readout.table <- rbindlist(list(readout.table, joined))
  }
  readout.table
}

parse.uniform <- function(readout.table)
{
  setDT(readout.table)[, Control := 0]
  setDT(readout.table)[Control == 0, Info := "regsample"]

  setDT(readout.table)[startsWith(GeneSymbol, "NTP"), GeneSymbol := "Scrambled"]
  setDT(readout.table)[GeneSymbol == "Scrambled", Control := -1]
  setDT(readout.table)[GeneSymbol == "Scrambled", Info:="negcontrol"]

  setDT(readout.table)[startsWith(GeneSymbol, "SARS"), GeneSymbol := gsub(" ", "",  GeneSymbol)]
  setDT(readout.table)[startsWith(GeneSymbol, "SARS"), Control := 1]
  setDT(readout.table)[startsWith(GeneSymbol, "SARS"), Info:="poscontrol"]

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

# get raw data for verification,
read.raw.data <- function(readout.data)
{
  files <- list.files(readout.data, full.names=T)
  full.raw.sars.table <- data.table()
  for (i in  seq(length(files)))
  {
    xls <- files[i]
    raw.sars.data <- as.data.table(read.xls(xls,sheet=7))
    raw.sars.data[,
                  which(
                    unlist(
                      lapply(raw.sars.data,
                             function(x) all(is.na(x)))))
                  :=NULL,
                  with=F]
    colnames(raw.sars.data) <- c("Plate",
                                 rep("Triplicate1", 12),
                                 rep("Triplicate2", 12),
                                 rep("Triplicate3", 12))
    plates <- as.integer(gsub("[^0-9]" ,"", raw.sars.data$Plate))
    ind <- which(!is.na(plates))
    plates <- rep(plates[ind], times=diff(c(ind, length(plates)+1)))
    raw.sars.data[,Plate := plates]
    formatted.sars.data <- raw.sars.data[,.SD,.SDcols=1:13]
    a <- cbind(1, raw.sars.data[,.SD,.SDcols=1:13])
    b <- cbind(2, raw.sars.data[,.SD,.SDcols=c(1,14:25)])
    c <- cbind(3, raw.sars.data[,.SD,.SDcols=c(1,26:37)])
    colnames(a) <- colnames(b) <- colnames(c) <-  c("Triplicate", "Plate", 1:12)
    formatted.sars.data <-  rbind(a,b,c)
    formatted.sars.data <- formatted.sars.data %>%
      group_by(Triplicate, Plate) %>% mutate(RowIdx=1:8)

    formatted.sars.data <- as.data.frame(dplyr::ungroup(formatted.sars.data))
    final.data <- data.table()
    n.row <- nrow(formatted.sars.data)
    n.col <- ncol(formatted.sars.data)
    for (ridx in seq(n.row))
    {
      rep <- formatted.sars.data$Triplicate[ridx]
      plate <- formatted.sars.data$Plate[ridx]
      rowidx <- formatted.sars.data$RowIdx[ridx]
      for (col in 3:14)
      {
        if (ridx == 1 & col == 3){
          final.data <- as.data.table(list(Replicate=i, Triplicate=rep,
                                           Plate=plate, RowIdx=rowidx,
                                           ColIdx=col-2,
                                           Readout=formatted.sars.data[ridx,col]))
        }
        else {
          final.data <- rbindlist(list(final.data,
                                       list(Replicate=i, rep, plate,
                                            rowidx, col-2,
                                            formatted.sars.data[ridx,col])))
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
      full.raw.sars.table <- final.data
    } else {
      full.raw.sars.table <- rbindlist(list(full.raw.sars.table, final.data))
    }
  }

  full.raw.sars.table
}

add.raw.data.to.readout <- function(readout.table, raw.data.table)
{
  n.row <- nrow(raw.data.table)
  n.row.prev <- nrow(readout.table)
  amb.wells <- c()
  sucs <- fls <- 0
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
                                   repl, as.character(NA), as.character(NA))))
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

# read readout
read.primary.readout <- function(readout.data)
{
  readout.table  <- get.readout.table(readout.data)
  readout.table <- parse.uniform(readout.table)
  raw.readout.table <- read.raw.data(readout.data)
  readout.table <- add.raw.data.to.readout(readout.table, raw.readout.table)
  readout.table
}

# read viability
read.viability <- function(viability.data)
{
  viability.table  <- get.readout.table(viability.data)
  viability.table  <- filter(viability.table, !is.na(Plate) & !is.na(Well)) %>%
                       dplyr::select(Replicate, Plate, Well, GeneSymbol,
                          RowIdx, ColIdx, Tri1, Tri2, Tri3)
  viability.table <- parse.uniform(viability.table)
  viability.table <- dplyr::select(viability.table, Replicate, Plate, Well,
                            GeneSymbol, RowIdx, ColIdx,
                            Tri1, Tri2, Tri3)
  raw.readout.table <- read.raw.data(viability.data)
  viability.table <- add.raw.data.to.readout(viability.table, raw.readout.table)
  colnames(viability.table)[which(startsWith(colnames(viability.table), "Tri"))] <-
    paste("Viability", 1:3, sep="")
  viability.table
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

  setDT(merged.layout.table)[,GeneSymbol := as.character(GeneSymbol)]
  setDT(merged.layout.table)[,PoolID := as.character(PoolID)]
  setDT(merged.layout.table)[,RefSeq := as.character(RefSeq)]

  setDT(readout.table)[,RefSeq := as.character(RefSeq)]
  setDT(readout.table)[,GeneSymbol := as.character(GeneSymbol)]
  setDT(readout.table)[,PoolID := as.character(PoolID)]
  setDT(readout.table)[,LocusID := as.character(LocusID)]

  final.table <- as.data.table(full_join(readout.table,
                  merged.layout.table,
                  by=c("RowIdx", "ColIdx", "PoolID")))

  setDT(final.table)[GeneSymbol.x == GeneSymbol.y,  GeneSymbol.y:=as.character(NA)]
  setDT(final.table)[RefSeq.x == RefSeq.y,  RefSeq.y:=as.character(NA)]
  final.table <- filter(final.table, !is.na(Plate) & !is.na(Well))
  colnames(final.table)[which(startsWith(colnames(final.table), "Tri"))] <-
    paste("Readout", 1:3, sep="")
  colnames(final.table)[which(startsWith(colnames(final.table), "Gene"))] <-
    c("GeneSymbol","AltGeneSymbol")
  colnames(final.table)[which(startsWith(colnames(final.table), "RefSeq"))] <-
    c("RefSeq","AltRefSeq")
  final.table
}

set.gene.name <- function(gene.x, gene.y)
{
  genes <- rep(NA_character_, length(gene.x))
  for (i in seq(length(gene.x)))
  {
    if (gene.x[i] == gene.y[i]) {gene <- gene.x[i]}
    else if (gene.x[i] %in% c("undef", "empty", NA)) {
      warning(paste("X", gene.x[i], "Y", gene[i],". Taking y"))
      gene <- gene.y[i]
    }
    else if (gene.y[i] %in% c("undef", "empty", NA))
    {
      warning(paste("X", gene.x[i], "Y", gene[i],". Taking x"))
      gene <- gene.x[i]
    }
    else
    {
      warning(paste("Ambiguous genes ", gene.x[i], gene.y[i]))
      gene <- paste(gene.x[i], gene.y[i], sep="_")
    }
    genes[i] <- gene
  }
 genes
}

join.readout.viability <-
function(readout.table, viability.table)
{
  final.table <- as.data.table(full_join(readout.table,
                                viability.table,
                                by=c("Replicate", "Plate", "Well", "RowIdx",
                                     "ColIdx")))
  setDT(final.table)[startsWith(GeneSymbol.x, "buffer"), GeneSymbol.x:="buffer"]
  setDT(final.table)[startsWith(GeneSymbol.y, "buffer"), GeneSymbol.y:="buffer"]
  final.table <- dplyr::mutate(final.table, GeneSymbol=set.gene.name(GeneSymbol.x, GeneSymbol.y))
  setDT(final.table)[,GeneSymbol.x := NULL]
  setDT(final.table)[, GeneSymbol.y := NULL]
  final.table[, 'Library':= 'DharmaconSMARTPool']
  setcolorder(final.table,
              c("Replicate", "Well", "RowIdx", "ColIdx", "Plate", "GeneSymbol",
                paste("Readout",1:3,sep=""), paste("Viability",1:3,sep=""),
                "Control", "Info", "Library", "AltGeneSymbol", "AltRefSeq",
                "LocusID", "RefSeq", "PoolID",
                "Sequences", "siRNAIDs"))
  final.table <- final.table[order(Well, Plate, GeneSymbol)]
  final.table
}

post.process <- function(screen)
{
  setDT(screen)[, Control := as.numeric(Control)]

  setDT(screen)[ColIdx==1, GeneSymbol := "undef"]

  setDT(screen)[RowIdx == 1 & ColIdx==1& Replicate==3, GeneSymbol := "undef"]
  setDT(screen)[RowIdx == 2 & ColIdx==1& Replicate==3, GeneSymbol := "buffer"]
  setDT(screen)[RowIdx == 3 & ColIdx==1& Replicate==3, GeneSymbol := "Scrambled"]
  setDT(screen)[RowIdx == 4 & ColIdx==1& Replicate==3, GeneSymbol := "buffer"]
  setDT(screen)[RowIdx == 5 & ColIdx==1& Replicate==3, GeneSymbol := "Scrambled"]
  setDT(screen)[RowIdx == 6 & ColIdx==1& Replicate==3, GeneSymbol := "buffer"]
  setDT(screen)[RowIdx == 7 & ColIdx==1& Replicate==3, GeneSymbol := "Scrambled"]
  setDT(screen)[RowIdx == 8 & ColIdx==1& Replicate==3, GeneSymbol := "EAV-nsp7"]

  setDT(screen)[(Plate == 1|Plate == 6) & RowIdx ==1 & ColIdx==1 & (Replicate==1|Replicate==2), GeneSymbol := "undef"]
  setDT(screen)[(Plate == 1|Plate == 6) & RowIdx ==2 & ColIdx==1 & (Replicate==1|Replicate==2), GeneSymbol := "buffer"]
  setDT(screen)[(Plate == 1|Plate == 6) & RowIdx ==3 & ColIdx==1 & (Replicate==1|Replicate==2), GeneSymbol := "Scrambled"]
  setDT(screen)[(Plate == 1|Plate == 6) & RowIdx ==4 & ColIdx==1 & (Replicate==1|Replicate==2), GeneSymbol := "GAPDH"]
  setDT(screen)[(Plate == 1|Plate == 6) & RowIdx ==5 & ColIdx==1 & (Replicate==1|Replicate==2), GeneSymbol := "SARS-nsp12 / 1"]
  setDT(screen)[(Plate == 1|Plate == 6) & RowIdx ==6 & ColIdx==1 & (Replicate==1|Replicate==2), GeneSymbol := "SARS-nsp12 / 2"]
  setDT(screen)[(Plate == 1|Plate == 6) & RowIdx ==7 & ColIdx==1 & (Replicate==1|Replicate==2), GeneSymbol := "SARS-nsp12 / 1 + 2"]
  setDT(screen)[(Plate == 1|Plate == 6) & RowIdx ==8 & ColIdx==1 & (Replicate==1|Replicate==2), GeneSymbol := "buffer"]

  setDT(screen)[startsWith(GeneSymbol, "NTP"), GeneSymbol := "Scrambled"]
  setDT(screen)[GeneSymbol == "Scrambled", Control := -1]
  setDT(screen)[GeneSymbol == "Scrambled", Info:="negcontrol"]

  setDT(screen)[startsWith(GeneSymbol, "SARS"), GeneSymbol := gsub(" ", "",  GeneSymbol)]
  setDT(screen)[startsWith(GeneSymbol, "SARS"), Control := 1]
  setDT(screen)[startsWith(GeneSymbol, "SARS"), Info:="poscontrol"]

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

layout.table                 <- read.sars.layout(lib.data)
primary.readout.table        <- read.primary.readout(readout.data)
joined.primary.layout.table  <-
  join.readout.layout(primary.readout.table,
                      layout.table)
viability.table <- read.viability(viability.data)
final.table <- join.readout.viability(joined.primary.layout.table,
                                      viability.table)
sars.kinase.screen <- final.table
sars.kinase.screen <- post.process(sars.kinase.screen)


save(file="inst/extdata/sars.kinase.screen.rda", sars.kinase.screen, compress="bzip2")
  write.table(
    sars.kinase.screen,
    file="~/PHD/data/data/sysvirdrug/integrated_data_files/sars_kinase_screen.tsv",
    sep='\t', quote=F, row.names=F)

#delete all and load data
 rm(list=ls())
 load('inst/extdata/sars.kinase.screen.rda')
