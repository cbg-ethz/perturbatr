# knockdown: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2015 - 2016 Simon Dirmeier
#
# This file is part of knockdown
#
# knockdown is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# knockdown is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with knockdown. If not, see <http://www.gnu.org/licenses/>.


#' @include class_knockdown_data.R
#' @include class_knockdown_analysed.R

#' @title Calculate statistics based on the hypergeometric-distribution to
#'  analyse the data.
#'
#' @description TODO
#'
#' @export
#' @docType methods
#' @rdname hyper_statistic-methods
#'
#' @param obj  the data to be analysed
#' @param padjust  multiple testing correction method
#' @param summ.method  summarize single siRNAs using mean or median
#' @param level  do hypergeometric test on gene level or siRNA level.
#'  If level=sirna is chosen, multiple replicates have to be given.
#'  If level=gene we use all the sirnas for a gene and optionally summarize
#'  siRNAs over replicate level.
#' @param do.summarization  boolean flag whether sirnas should be summarized
#'  if level=gene is chosen
#' @param hit.ratio  ratio of succesful hits such that a gene is considered a
#'  'hit'
#' @param effect.size  the mean effect size used for hit prioritization
#' @param pval.threshold  p-value threshold when a siRNA should be considered a
#'  'hit'
#' @param qval.threshold  q-value threshold when a siRNA should be considered a
#'  'hit'
#'
#' @return returns a \code{knockdown.hyper.analysed} object
#'
setGeneric(
  "hyper.statistic",
  function(obj,
           padjust          = c("BH", "bonferroni"),
           summ.method      = c("mean", "median"),
           level            = c("gene", "sirna"),
           do.summarization = FALSE,
           hit.ratio        = 0.5,
           effect.size      = 0,
           pval.threshold   = 0.05,
           qval.threshold   = 1)
  {
    standardGeneric("hyper.statistic")
  },
  package="knockdown"
)

#' @rdname hyper_statistic-methods
#' @aliases hyper.statistic,knockdown.data-method
#' @import data.table
#' @importFrom methods new
setMethod(
  "hyper.statistic",
  signature = signature(obj="knockdown.normalized.data"),
  function(obj,
           padjust          = c("BH", "bonferroni"),
           summ.method      = c("mean", "median"),
           level            = c("gene", "sirna"),
           do.summarization = FALSE,
           hit.ratio        = 0.5,
           effect.size      = 0,
           pval.threshold   = 0.05,
           qval.threshold   = 1)
  {
    stopifnot(is.logical(do.summarization))

    res <- .hyper.statistic(
      obj              = obj@.data,
      padjust          = match.arg(padjust),
      summ.method      = match.arg(summ.method),
      level            = match.arg(level),
      do.summarization = do.summarization)
    priorit <- .prioritize.hyper.statistic(
      res, hit.ratio, effect.size, pval.threshold, qval.threshold)

    ret <- methods::new(
      "knockdown.hyper.analysed",
      .gene.hits = data.table::as.data.table(priorit),
      .data      = data.table::as.data.table(obj@.data),
      .params=list(effect.size=effect.size,
                   hit.ratio=hit.ratio,
                   pval.threshold=pval.threshold,
                   qval.threshold=qval.threshold))
    ret
  }
)

#' @noRd
#' @import data.table
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
.hyper.statistic <- function(obj,
                             padjust,
                             summ.method,
                             level,
                             do.summarization)
{
  if (do.summarization & level=="sirna")
    stop("Cant do summarization on sirna level. Choose level=gene")
  if (level=="gene" & do.summarization)
    message(paste("Summarizing with ", summ.method, "!", sep=""))

  summ.method <- .summarization.method(summ.method)
  ret <- dplyr::group_by(obj, Virus, Screen, ReadoutType,
                         ScreenType, Library, Design, Cell) %>%
    dplyr::mutate(grp = .GRP) %>% ungroup
  grps <- unique(ret$grp)

  res <- data.table::rbindlist(
    lapply(
      grps, function (g)
      {
        grp.dat <- dplyr::filter(ret, grp==g)
        message(paste("Doing grp: ", paste(grp.dat$Virus[1],
                                           grp.dat$Screen[1],
                                           grp.dat$ScreenType[1],
                                           grp.dat$ReadoutType[1],
                                           grp.dat$Cell[1],
                                           grp.dat$Design[1],
                                           grp.dat$Library[1],
                                           sep=", ")))
        fr <- .do.hyperstatistic(grp.dat, padjust, summ.method,
                                 do.summarization, level)
        fr
      }
    )
  )

 res
}

#' @noRd
#' @import data.table
#' @importFrom dplyr select filter group_by summarize mutate
#' @importFrom tidyr separate
#' @importFrom stats p.adjust
.do.hyperstatistic <- function(obj, padjust, summ.method,
                               do.summarization, level)
{
  res <- obj %>% ungroup
  if (obj$Design[1] == "single" & level=="gene" & do.summarization)
  {
    message(paste("\t...summarizing single siRNAs over replicates!"))
    # summarize all the sirnas over the different replicates
    # so: for a gene A and siRNA B, every siRNA B has ONE observation
    res <- dplyr::group_by(res, Virus, Screen, Library,
                           ScreenType, ReadoutType,
                           Cell, Design,
                           Plate, RowIdx, ColIdx,
                           GeneSymbol, Entrez, siRNAIDs) %>%
      dplyr::summarize(Readout=summ.method(Readout, na.rm=TRUE)) %>%
      ungroup
  }
  else
  {
    message(paste("\t...NOT summarizing single siRNAs over replicates!"))
  }
  if (obj$Design[2] == "pooled")
  {
    level <- "gene"
    message("\t...setting level=gene since a pooled library is used!")
  }

  # do hyper test on every screen
  res <- dplyr::group_by(res, Virus, Screen, Library, ReadoutType, ScreenType,
                         Cell, Design) %>%
    dplyr::mutate(HRes=.hypertest(GeneSymbol, siRNAIDs, Plate,
                                  RowIdx, ColIdx, Readout, level)) %>%
    ungroup %>%
    tidyr::separate(HRes, c("Pval", "Hit"), sep="_")

  if ("grp" %in% colnames(res)) res <- dplyr::select(res, -grp)

  data.table::setDT(res)[,Pval := as.numeric(Pval)]
  data.table::setDT(res)[,Qval := p.adjust(Pval, method=padjust)]
  data.table::setDT(res)[,Hit  := as.logical(as.numeric(Hit))]
  res <- res[order(Pval)]
  data.table::setDT(res)[,HyperRank := cumsum(res$Hit)]
  data.table::setDT(res)[Hit == 0, HyperRank := NA_integer_]
  assertthat::assert_that(all(order(res$Pval) == order(res$Qval)))

  res
}

#' @noRd
#' @import data.table
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
.hypertest <- function(all.genes, sirnas, plates, rows, cols, readouts, level)
{
  fr <- data.table::data.table(genes=all.genes, sirnas=sirnas,
                               plates=plates, rows=rows, cols=cols,
                               # leave the abs, since we rank the whole thing
                               readout=abs(readouts),
                               ord=1:length(all.genes)) %>%
    dplyr::mutate(rank=rank(-readout, ties.method="max")) %>%
    .[order(rank)]

  ## this part is tricky!
  ## if we do the hypergeometric test on genes we need another grouping
  ## as for siRNAs
  # on gene level group by genes and take all siRNAs for test
  if (level == "gene")
  {
    fr <- dplyr::group_by(fr, genes)
  }
  # on sirna level group by gene name and specific sirna
  # (the siRNA SHOULD be identical to grouping by plate/row/col;
  #  so the last grouping should be redundant)
  else if (level == "sirna")
  {
    fr <- dplyr::group_by(fr, genes, sirnas, plates, rows, cols)
  }
  else
  {
    stop("Please provide a standard method!")
  }

  # used the grouped ranks and all genes and to hyper test
  fr <- fr %>%
    dplyr::mutate(H=.hypertest.grp(rank, all.genes)) %>% ungroup %>%
    .[order(ord)]

  fr$H
}

#' Test as described in the original RSA paper by Koenig et al, 2007
#' @noRd
.hypertest.grp <- function(ranks, all.genes)
{
  N      <- length(all.genes)
  ndrawn <- length(ranks)

  hi       <- .hypertest.grp.test(ndrawn, ranks, N)
  min.idx  <- which.min(hi[, 1])
  min.prob <- hi[min.idx, 1]
  is.hit   <- as.numeric(seq(ndrawn) <= hi[min.idx, 2])

  paste(min.prob, is.hit, sep="_")
}

#' @noRd
#' @importFrom stats phyper
.hypertest.grp.test <- function(ndrawn, ranks, N)
{
  t(sapply(1:ndrawn, function(i)
  {
    prob <- stats::phyper(i - 1, ranks[i], N - ranks[i],
                          ndrawn, lower.tail = FALSE, log.p=FALSE)
    prob <- max(prob, 0.0)
    cutoff <- i
    c(prob=prob, cutoff=cutoff)
  }))

}

#' @noRd
#' @import data.table
#' @importFrom metap sumlog
#' @importFrom dplyr group_by summarize ungroup filter select
.prioritize.hyper.statistic <- function(obj,
                                        hit.ratio,
                                        effect.size,
                                        pval.threshold,
                                        qval.threshold)
{
  if (effect.size != 0)    stop("Effect size not yet implemented")
  if (qval.threshold != 1) stop("Q-value not  yet implemented")
  res <- dplyr::group_by(obj, Virus, Screen, Library,
                         ScreenType, ReadoutType,
                         Design, Cell,
                         GeneSymbol, Entrez) %>%
    dplyr::summarize(HitRatio   = (sum(Hit == TRUE, na.rm=TRUE) / n()),
                     Pval  = metap::sumlog(Pval)$p,
                     Qval  = metap::sumlog(Qval)$p,
                     MeanEffect = mean(Readout,na.rm=TRUE),
                     MaxEffect  = max(Readout, na.rm=TRUE),
                     MinEffect  = min(Readout, na.rm=TRUE),
                     MinPval    = min(Pval, na.rm=TRUE),
                     MinQval    = min(Qval, na.rm=TRUE),
                     AllPval=paste(sprintf("%03f", Pval), collapse=","),
                     AllQval=paste(sprintf("%03f", Qval), collapse=",")) %>%
    ungroup %>%
    dplyr::filter(HitRatio >= hit.ratio)

  res
}
