norm <- function(read, rows, cols, ge, cont)
{
  r <- read - mean(read[cols == max(cols)])
  (r / mean(r[which(cont == -1 & ge == "Scrambled")]) * 100)
}


tet <- function(read, gen)
{
  idx <- which(gen == "Scrambled")
  #maybe what is done here
  (read)
}


test.st <- function(read, mu)
{
  t.test(read, mu=mu)$p.value
}

# dat is the raw CHIKV dataset
stest <-
  filter(dat, ReadoutClass=="Readout", Replicate==1, Idx==3) %>%
  group_by(Replicate, Plate, Idx) %>%
  dplyr::mutate(Readout = norm(Readout, RowIdx, ColIdx, GeneSymbol, Control))

# first normalize:
#   substract mean of last column
#   then divide by mean of Scrambled DNA
# second summarize triplicates using mean
# third summarize replicates using mean
# then mutate every gene using log fold change
# finally calculate t statistics using
s <-
  filter(dat, ReadoutClass=="Readout") %>%
  group_by(Replicate, Plate, Idx) %>%
  dplyr::mutate(Readout = norm(Readout, RowIdx, ColIdx, GeneSymbol, Control)) %>%
  ungroup %>%
  dplyr::group_by(Replicate, ReadoutClass, GeneSymbol, Plate, RowIdx, ColIdx) %>%
  dplyr::summarize(Readout=mean(Readout)) %>% .[order(GeneSymbol)] %>%
  ungroup %>%
  group_by(GeneSymbol) %>%
  dplyr::summarize(Readout=mean(Readout)) %>% ungroup %>%
  dplyr::mutate(Readout=tet(Readout, GeneSymbol)) %>% .[order(Readout, decreasing=F)]

#mean value for t-test (scrambled rna means)
mu <-
  filter(dat, ReadoutClass=="Readout") %>%
  group_by(Replicate, Plate, Idx) %>%
  dplyr::mutate(Readout = norm(Readout, RowIdx, ColIdx, GeneSymbol, Control)) %>%
  filter(GeneSymbol=="Scrambled") %>% ungroup %>%
  dplyr::summarize(Readout=mean(Readout)) %>% select(Readout) %>% unlist

# test significance of gene and correct using BH
test <-
  filter(dat, ReadoutClass=="Readout") %>%
  group_by(Replicate, Plate, Idx) %>%
  dplyr::mutate(Readout = norm(Readout, RowIdx, ColIdx, GeneSymbol, Control)) %>%
  dplyr::group_by(Replicate, ReadoutClass, GeneSymbol, Plate, RowIdx, ColIdx) %>%
  dplyr::summarize(Readout=mean(Readout)) %>% .[order(GeneSymbol)] %>%
  dplyr::group_by(GeneSymbol) %>%
  dplyr::summarize(pval=test.st(Readout, mu)) %>% .[order(pval)] %>% ungroup %>%
  dplyr::mutate(pval=p.adjust(pval, "BH"))

res <- full_join(s, test) %>% filter(pval < .05)
