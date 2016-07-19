#' Load the pre-calculated hugo2Entrez mapping
#'
#' @noRd
#' @importFrom data.table setDT
load.hugo.mappings <-
  function ()
  {
    hugo2entrez <-
      as.data.table(read.table("~/PHD/data/data/sysvirdrug/mappings/geneIDmappings/hugo2entrez_screen_final.tsv",
                               sep="\t", fill=T, header=F))
    colnames(hugo2entrez) <- c("GeneSymbol", "Entrez", "Source")
    setDT(hugo2entrez)[ ,Source := NULL]
    hugo2entrez$GeneSymbol <- as.character(hugo2entrez$GeneSymbol)
    hugo2entrez
  }

#' Correct some nasty gene names
#'
#' @noRd
#' @importFrom data.table setDT
annotate <- function(rnai.screen.summarized)
{
  setDT(rnai.screen.summarized)[GeneSymbol == "Trim5alpha", GeneSymbol := "TRIM5"]
  setDT(rnai.screen.summarized)[GeneSymbol == "p53", GeneSymbol := "TP53"]
  setDT(rnai.screen.summarized)[GeneSymbol == "LEDGF-p75", GeneSymbol := "PSIP1"]
  setDT(rnai.screen.summarized)[GeneSymbol == "PP1A", GeneSymbol := "PPP1CA"]
  setDT(rnai.screen.summarized)[GeneSymbol == "PP1B", GeneSymbol := "PPP1CB"]
  setDT(rnai.screen.summarized)[GeneSymbol == "PP1C", GeneSymbol := "PPP1CC"]
  setDT(rnai.screen.summarized)[GeneSymbol == "CKI-delta", GeneSymbol := "CSNK1D"]
  setDT(rnai.screen.summarized)[GeneSymbol == "CKI-alpha", GeneSymbol := "CSNK1A1"]
  setDT(rnai.screen.summarized)[GeneSymbol == "CKI-alpha", GeneSymbol := "CSNK1A1"]
  setDT(rnai.screen.summarized)[GeneSymbol == "bCop", GeneSymbol := "COPB"]
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
