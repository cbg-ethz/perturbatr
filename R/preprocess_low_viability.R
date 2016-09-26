# knockout: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2015 - 2016 Simon Dirmeier
#
# This file is part of knockout
#
# knockout is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# knockout is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with knockout. If not, see <http://www.gnu.org/licenses/>.

#' @noRd
#' @import data.table
#' @importFrom dplyr group_by mutate filter select
.rm.cytotoxic <-
function
(
  obj,
  rm.cytotoxic
)
{
  if (!is.null(rm.cytotoxic) & "Viability" %in% colnames(obj))
  {
    message(paste("Removing un-viable genes by comparing to ",
                  rm.cytotoxic, " viability.", sep=""))
    message(paste("\tUsing t-test with mu=85*mean(viability(",
                  rm.cytotoxic ,")) and alternative=less", sep=""))
    obj <-
      dplyr::group_by(obj, Virus, Screen, Library,
                      InfectionType, ReadoutType,
                      Design, Cell) %>%
      # set siRNAS that are cytotoxic within each replicate
      dplyr::mutate(Remove=.set.cytotoxic(Readout, Viability, siRNAIDs,
                                          Control, GeneSymbol,
                                          Plate, RowIdx, ColIdx,
                                          rm.cytotoxic)) %>%
      ungroup
  }
  invisible(obj)
}

#' Summarize the sirna readouts and then decide whether sirnas should be removed from the analysis.
#'
#' @noRd
#'
#' @import data.table
#' @importFrom assertthat assert_that
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr filter
.set.cytotoxic <-
function
(
  re,
  val,
  sirnas,
  ctrl,
  genes,
  plates,
  rows,
  cols,
  comp.to
)
{
  # against which gene should viability be compared
  comp.to <- tolower(comp.to)
  genes <- tolower(genes)
  remarr <- rep(F, length(re))
  if (!all(is.na(val)))
  {
    cont.idx <- which(ctrl == -1 & genes == comp.to)
    if (is.na(cont.idx[1])) stop("No control index found! Change gene names")
    assertthat::assert_that(length(cont.idx) == length(which(genes == comp.to)))
    # get mean readout and  viability of scrambled siRNAs
    cont.vial.thresh <- base::mean(val[cont.idx], na.rm=T) * .85
    cont.re   <- base::mean(re [cont.idx], na.rm=T)
    # get the mean readouts and viabilities for every siRNA
    fr <-
      data.table::data.table(Re=re, Val=val, Gene=genes, Sirna=sirnas,
                             Plate=plates, Row=rows, Col=cols, Control=ctrl) %>%
      # group every screen by gene, plate, sirna, row and column
      dplyr::group_by(Sirna, Gene, Plate, Row, Col, Control) %>%
      dplyr::mutate(n = n()) %>% ungroup
      if (any(fr$n < 3))
        message(paste("\tSome siRNAs have less than three observations.",
                      "Setting p.val=1 to prohibit removal of siRNAs erroneously."))
    fr <-
      dplyr::group_by(fr, Sirna, Gene, Plate, Row, Col, Control) %>%
      # summarize over replicates
      dplyr::summarize(MR=mean(Re, na.rm=T),
                       MV=mean(Val, na.rm=T),
                       Pval=.t.test.vial(Val, cont.vial.thresh, Gene,
                                         Sirna, Plate, Row, Col)) %>%
      ungroup
    # this .85 is magic
    # get the sirnas that show toxicity
    lowe <- dplyr::filter(fr,
                          MV < cont.vial.thresh,
                          MR < cont.re,
                          Pval < .05,
                          !is.na(Sirna),
                          !is.na(Gene),
                          Control == 0,
                          Gene != "scrambled",
                          Gene != "buffer",
                          Gene != "empty",
                          Gene != "gapdh",
                          Gene != "gfp")
    remarr[which(sirnas %in% lowe$Sirna &
                   rows %in% lowe$Row&
                   plates %in% lowe$Plate &
                   cols %in% lowe$Col)
           ] <- T
  }
  remarr
}

#' @noRd
#' @importFrom stats t.test
.t.test.vial <-
function
(
  vial,
  cont.vial.thresh,
  gene, sirna, plate,
  row, col
)
{
  p.val <- 1.0
  if (is.na(gene) | is.na(sirna)) p.val <- 1
  else if (length(p.val) < 3)     p.val <- 1
  else if (!all(is.na(vial)))
  {
    war <- paste("Plate:", plate, ", row:", row, ", col:", col,
                 ", gene:", gene, ", sirna:", sirna, "!", sep="" )
    tryCatch({
      p.val <- t.test(vial, mu=cont.vial.thresh, alternative="less")$p.value },
      warning = function(r) { message(paste(r, war)) },
      error   = function(r) { message(paste(r, war)) }
    )
  }
  p.val
}
