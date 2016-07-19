#################################
#
# read DENV kinase data
# this is a CONST FILE reader and should actually not be used since the data is also in data/denv.kinase.screen.Rd

require(data.table)
require(gdata)
rm(list=ls())

den.dir <- '/Users/simondi/PHD/data/data/sysvirdrug/screening_data/dengue_kinase_screen/'
lib.data <- paste(den.dir,  'ambion_kinase_library_layout/', sep='')
readout.data <- paste(den.dir,  'DENV_kinase_screen_results/', sep='')

# read mappings
files <- list.files(lib.data)
den.kinase.ambion.library.layout <- list()
for (file in files)
{
  idx <- as.numeric(gsub("[^0-9]", "", file))
  mapping <- read.csv(paste(lib.data, file, sep=''), sep='\t')
  colnames(mapping) <- c("WellIdx", "Row", "Col", "NoIdea", "siRNA", "GeneSymbol")
  den.kinase.ambion.library.layout[[idx]] <- mapping
}

# read siRNA info
den.kinase.siRNA.info <- read.csv(paste(den.dir, "kinases_ambion_siRNA_info.txt", sep=''), sep='\t')[, c(1,2,3,5,6,8)]

# read readout
files <- list.files(readout.data, full.names=T)
files <-  files[grep("*/e\\d+l\\d+.txt", files)]
readout.table <- matrix()
for (i in  seq(length(files)))
{
  f <- files[i]
  basename <-  sub('.txt', '',  basename(f))
  ve <- na.omit(as.numeric(unlist(base::strsplit(gsub("[^0-9]", ' ', basename), " "))))
  readout <- read.csv(f, sep='\t')[,c(1,3,9)]
  if (i == 1) readout.table <- cbind(matrix(ve, nrow(readout), length(ve),byrow=T), readout)
  else
  {
    mat <- cbind(matrix(ve, nrow(readout), length(ve), byrow=T), readout)
    readout.table <- rbind(readout.table, mat)
  }
}
colnames(readout.table)[1:2] <- c("Replicate", "Plate")
readout.table <- as.data.table(readout.table[order(readout.table$Replicate, readout.table$Plate), ])

# add sirna and plate meta data to every row
for (i in seq(length(den.kinase.ambion.library.layout)))
{
  plate <- den.kinase.ambion.library.layout[[i]]
  for (j in seq(dim(plate)[1]))
  {
    plate.meta.info <- plate[j, c(2,3,5,6)]
    well.idxs <- which(readout.table$Plate==i & readout.table$WellNo==j)
    readout.table[well.idxs, c("RowIdx", "ColIdx", "siRNA", "GeneSymbol"):=plate.meta.info]
  }
}

# add gene RefSeq and siRNA sequence
for (i in seq(length(den.kinase.siRNA.info[,1])))
{
  sirna.info <- den.kinase.siRNA.info[i, ]
  idx <- which(readout.table$siRNA == unlist(sirna.info[1]))
  readout.table[idx, c('RefSeq', 'Sequence') := sirna.info[c(4,6)] ]
}

readout.table[, 'Library' := 'Ambion']
readout.table[, 'Control' := 0]
# add information about controls
setDT(readout.table)[Control == 0, Info := "regsample"]
setDT(readout.table)[GeneSymbol == "Scrambled", Control := -1]
setDT(readout.table)[GeneSymbol == "Scrambled", Info:="negcontrol"]

setDT(readout.table)[GeneSymbol %in% c('DV-E', 'DV-NS3', "DV-NS5"), Control := 1]
setDT(readout.table)[GeneSymbol %in% c('DV-E', 'DV-NS3', "DV-NS5"), Info:="poscontrol"]

setDT(readout.table)[startsWith(GeneSymbol, "HCV"), Control := -1]
setDT(readout.table)[startsWith(GeneSymbol, "HCV"), Info := "negcontrol"]

setDT(readout.table)[startsWith(GeneSymbol, "empty"), Control := NA_integer_]
setDT(readout.table)[startsWith(GeneSymbol, "empty"), Info := NA_character_]
setDT(readout.table)[startsWith(GeneSymbol, "empty"), GeneSymbol := NA_character_]
setDT(readout.table)[startsWith(siRNA, "empty"), siRNA := NA_character_]

# set missing values
setDT(readout.table)[RefSeq == '', RefSeq := as.character(NA)]

denv.kinase.screen <- readout.table
save(file="inst/extdata/denv.kinase.screen.rda", denv.kinase.screen, compress="bzip2")
write.table(denv.kinase.screen, file="~/PHD/data/data/sysvirdrug/integrated_data_files/denv_kinase_screen.tsv",sep='\t',quote=F,row.names=F)

# delete all and load data
rm(list=ls())
load('inst/extdata/denv.kinase.screen.rda')
