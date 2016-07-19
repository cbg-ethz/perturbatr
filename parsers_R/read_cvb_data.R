#################################
#
# read CVB3 data
# this is a CONST FILE reader and should actually not be used since the data is also in data/cvb.replicon.screen.Rd

require(data.table)
require(dplyr)
require(reshape2)
require(tidyr)
require(gdata)
library(R.matlab)
rm(list=ls())


cvb.dir      <- '/Users/simondi/PHD/data/data/sysvirdrug/screening_data/cvb_screen/'
readout.file <- paste(cvb.dir,  'cvb_full_data.mat', sep='')

matlab.content <- R.matlab::readMat(readout.file)
cvb.data <- matlab.content$BASICDATA
names <- rownames(cvb.data)
name.map <- list()
for (i in seq(length(names)))
{
  name.map[[names[i]]] <- i
}
# 75 x 384 matrix ... so 75 plates!?
cell.counts <- cvb.data[[1]]


## TODO THIS IS NOT DONE YET BECAUSE DATA IS MISSING
cvb.replicon.screen <- readout.table
save(file="data/cvb.replicon.screen.rda", cvb.replicon.screen, compress="bzip2")
write.table(cvb.replicon.screen,
            file="~/PHD/data/data/sysvirdrug/integrated_data_files/cvb_replicon_screen.tsv",
            sep='\t',quote=F,row.names=F)

#delete all and load data
rm(list=ls())
load('data/cvb.replicon.screen.rda')

