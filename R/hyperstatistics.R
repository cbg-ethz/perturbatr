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

#' Calculate statistics based on the hypergeometric-distribution to analyse the data.
#'
#' For this you should use a standardization before.
#'
#' @export
#' @import data.table
#'
#' @param obj  the data to be analysed
#' @param padjust  multiple testing correction method
#' @param summ.method  summarize single siRNAs using mean or median
#' @param level  do hypergeometric test on gene level or siRNA level
#' @param ...   additional params
hyperstatistic <-
function
(
  obj,
  padjust=c("BH", "bonferroni"),
  summ.method=c("mean", "median"),
  level=c("gene", "sirna"),
  ...
)
{
  UseMethod("hyperstatistic", obj)
}

#' @noRd
#' @export
#' @import data.table
hyperstatistic.svd.data <-
function
(
  obj,
  padjust=c("BH", "bonferroni"),
  summ.method=c("mean", "median"),
  level=c("gene", "sirna"),
  ...
)
{
  summ.method <- match.arg(summ.method)
  padjust     <- match.arg(padjust)
  level       <- match.arg(level)
  ret <- .hyperstatistic(obj, padjust=padjust,
                         summ.method=summ.method, level=level, ...)
  class(ret) <- c("svd.analysed.hyper", "svd.analysed", class(ret))
  invisible(ret)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
.hyperstatistic <-
function
(
  obj,
  padjust,
  summ.method,
  level,
  ...
)
{
  message(paste("Correcting with ", padjust, "!", sep=""))
  if (level=="gene" & !is.na(summ.method))
    message(paste("Summarizing with ", summ.method, "!", sep=""))
  summ.method <- .summarization.method(summ.method)
  grp.indexes <- dplyr::group_indices(obj, Virus, Screen, ReadoutType,
                                      InfectionType, Library, Design, Cell)
  ret <-  dplyr::mutate(obj, grp=grp.indexes)
  grps <- unique(ret$grp)
  res <- do.call(
    "rbind",
    lapply
    (
      grps,
      function (g)
      {
        grp.dat <- dplyr::filter(ret, grp==g)
        message(paste("Doing grp: ", paste(grp.dat$Virus[1],
                                         grp.dat$Screen[1],
                                         grp.dat$InfectionType[1],
                                         grp.dat$ReadoutType[1],
                                         grp.dat$Cell[1],
                                         grp.dat$Design[1],
                                         grp.dat$Library[1],
                                         sep=", ")))
        fr <- .do.hyperstatistic(grp.dat, padjust, summ.method, level)
        fr
      }
    )
  )
  invisible(res)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr mutate
#' @importFrom tidyr separate
#' @importFrom stats p.adjust
.do.hyperstatistic <-
function
(
  obj,
  padjust,
  summ.method,
  level
)
{
  # TODO: make this nicer and split up
  res <- obj %>% ungroup
  if (obj$Design[1] == "single" & level=="gene")
  {
    message(paste("\t..summarizing single siRNAs over replicates!"))
    res <- dplyr::group_by(res, Virus, Screen, Library,
                           InfectionType, ReadoutType,
                           Cell, Design,
                           Plate, RowIdx, ColIdx,
                           GeneSymbol, Entrez, siRNAIDs) %>%
      dplyr::summarize(Readout=summ.method(Readout, na.rm=T)) %>%
      ungroup
  } else {
    message(paste("\t..NOT summarizing single siRNAs over replicate!"))
  }
  if (obj$Design[2] == "pooled")
  {
    level <- "gene"
    message("\t...setting level=gene since a pooled library is used!")
  }
  res <-
    dplyr::group_by(res, Virus, Screen, Library,
                    ReadoutType, InfectionType, Cell, Design) %>%
    dplyr::mutate(HRes=.hypertest(GeneSymbol, siRNAIDs, Plate,
                                  RowIdx, ColIdx, Readout, level)) %>%
    ungroup %>%
    tidyr::separate(HRes, c("Pval", "Hit"), sep="_")
  if ("grp" %in% colnames(res)) res <- dplyr::select(res, -grp)
  data.table::setDT(res)[,Pval := as.numeric(Pval)]
  data.table::setDT(res)[,Hit := as.logical(as.numeric(Hit))]
  res <- res[order(Pval)]
  data.table::setDT(res)[,HyperRank := cumsum(res$Hit)]
  data.table::setDT(res)[Hit == 0, HyperRank := NA_integer_]
  res <- dplyr::mutate(res, Pvalcorr=p.adjust(Pval, method=padjust))
  invisible(res)
}

#' @noRd
#' @import data.table
#' @import dtplyr
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
.hypertest <-
function
(
 all.genes,
 sirnas,
 plates,
 rows,
 cols,
 readouts,
 level
)
{
  fr <- data.table::data.table(genes=all.genes, sirnas=sirnas,
                               plates=plates, rows=rows, cols=cols,
                               readout=abs(readouts),
                               ord=1:length(all.genes)) %>%
    dplyr::mutate(rank=rank(-readout, ties.method="max")) %>%
    .[order(rank)]
  if (level == "gene")
  {
    fr <- dplyr::group_by(fr, genes)
  }
  else if (level =="sirna")
  {
    fr <- dplyr::group_by(fr, genes, sirnas, plates, rows, cols)
  }
  else
  {
    stop("Please provide a standard method!")
  }
  fr <- fr %>%
    dplyr::mutate(H=.hypertest.grp(rank, all.genes)) %>%
    ungroup %>%
    .[order(ord)]
  fr$H
}

#' @noRd
#' @importFrom stats phyper
.hypertest.grp <-
function
(
  ranks,
  all.genes
)
{
  N <- length(all.genes)
  ndrawn <- length(ranks)
  hi <- t(sapply(1:ndrawn, function(i)
  {
    prob <- stats::phyper(i - 1, ranks[i], N - ranks[i],
                   ndrawn, lower.tail = F,log.p=F)
    prob <- max(prob, 0.0)
    cutoff <- i
    c(prob=prob, cutoff=cutoff)
  }))
  min.idx <- which.min(hi[,1])
  min.prob <- hi[min.idx,1]
  is.hit <- as.numeric(seq(ndrawn) <= hi[min.idx,2])
  paste(min.prob, is.hit, sep="_")
}
