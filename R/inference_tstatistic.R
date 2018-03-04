# perturbatr: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2018 Simon Dirmeier
#
# This file is part of perturbatr
#
# perturbatr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# perturbatr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with perturbatr. If not, see <http://www.gnu.org/licenses/>.


#' @include class_data.R
#' @include class_analysed.R


#' @title T-test
#'
#' @description First test statistics for perturbations
#'  based on Student's t distribution are computed. Then the test
#'  statistics for
#'  perturbations, like siRNAs or gRNAs, are aggregated to a signficance level
#'  for the genes that they correspond to. That means, for instance, for 10
#'  gRNAs that target gene A, first the statistical significance for the 10
#'  guides is estimated, and finally the significance of gene A is computed.
#'
#' @export
#' @docType methods
#' @rdname tStatistic-methods
#'
#' @import data.table
#'
#' @param obj  a \code{\link{PerturbationData}} object
#' @param mu  side to which the mean of a perturbation is compared to
#' @param padjust  multiple testing correction method
#' @param hit.ratio  the ratio of suffesfully identified perturbations for a
#'  gene, in order to make it count as statistically significant. For instance,
#'  for a CRISPR screen where 10 guides have been used for a gene, a
#'  \code{hit.ratio} of .5 would require 5 guide to be estimated significant.
#' @param effect.size  the relative strength of a signal to count as a hit,
#'  i. e. the biological significance of a gene/perturbation
#' @param pval.threshold  p-value threshold when a perturbation should be
#'  considered a 'hit', e.g. that the perturbation resulted in a significant
#'  change in phenotype
#' @param qval.threshold  q-value threshold when a perturbation should be
#'  considered a 'hit', e.g. that the perturbation resulted in a significant
#'  change in phenotype
#'
#' @return returns a \code{perturbation.tstatistic.analysed} object
#' @examples
#'  library(magrittr)
#'  data(rnaiscreen)
#'  v1.dat <- perturbatr::filter(rnaiscreen, Condition=="V1")
#'  v1.data.norm <- perturbatr::preprocess(v1.dat, normalize="z.score") %>%
#'    filter(GeneSymbol == "c1qb", Screen=="Kinome")
#'  v1.res <- tstatistic(v1.data.norm)
setGeneric(
  "tStatistic",
  function(obj,
           mu=c(0, "scrambled", "control"),
           padjust=c("BH", "bonferroni"),
           hit.ratio=0.5,
           effect.size=0,
           pval.threshold=0.05,
           qval.threshold=1)
  {
    standardGeneric("tstatistic")
  },
  package="perturbation"
)


#' @rdname t_statistic-methods
#' @aliases tstatistic,perturbation.data-method
#' @import data.table
#' @importFrom methods new
setMethod(
  "tstatistic",
  signature = signature(obj = "perturbation.normalized.data"),
  function(obj,
           mu=c(0, "scrambled", "control"),
           padjust=c("BH", "bonferroni"),
           hit.ratio=0.5,
           effect.size=0,
           pval.threshold=0.05,
           qval.threshold=1)
  {
    res <- .t.statistic(obj@.data,
                        mu      = match.arg(mu),
                        padjust = match.arg(padjust))
    priorit <- .prioritize.tstatistic(
      res, hit.ratio, effect.size, pval.threshold, qval.threshold)

    ret <- methods::new(
      "perturbation.tstatistic.analysed",
      .gene.hits = data.table::as.data.table(priorit),
      .data=res,
      .params=list(effect.size=effect.size,
                   hit.ratio=hit.ratio,
                   pval.threshold=pval.threshold,
                   qval.threshold=qval.threshold))
    ret
  }
)

#' @noRd
#' @import data.table
.t.statistic <- function(obj, mu, padjust)
{
  message(paste("Correcting with ", padjust, "!", sep=""))
  message(paste("Taking", mu, "for t-test mu!"))
  # do NOT group by replicate here, we do inference using these
  ret <- dplyr::group_by(obj, Condition, Screen, Library,
                         ReadoutType, ScreenType,
                         Design, Cell) %>%
    dplyr::mutate(grp=.GRP) %>%
  ungroup
  grps <- unique(ret$grp)
  res <- do.call(
  "rbind", lapply(
    grps, function (g)
    {
      grp.dat <- dplyr::filter(ret, grp==g) %>% ungroup
      message(paste("Doing grp: ", paste(grp.dat$Condition[1],
                                         grp.dat$Screen[1],
                                         grp.dat$ScreenType[1],
                                         grp.dat$ReadoutType[1],
                                         grp.dat$Cell[1],
                                         grp.dat$Design[1],
                                         grp.dat$Library[1],
                                         sep=", ")))
        fr <- .do.t.statistic(grp.dat, padjust, mu)
        fr
      }
    )
  )

  res
}

#' @noRd
#' @import data.table
#' @importFrom dplyr mutate group_by summarize ungroup
#' @importFrom stats p.adjust
#' @importFrom assertthat assert_that
.do.t.statistic <- function(obj, padjust, mu)
{
  res <- obj %>% ungroup %>%
    dplyr::group_by(Plate) %>% dplyr::mutate(grp=.GRP) %>% ungroup
  grps <- unique(res$grp)
  ret  <- do.call(
    "rbind",
    lapply(
      grps, function(g)
      {
        grp.dat <- dplyr::filter(res, grp==g) %>% ungroup
        fr      <- .t.statisic.plate(grp.dat, mu)
        fr
      }
    )
  )
  data.table::setDT(ret)[,Pval := as.numeric(Pval)]
  data.table::setDT(ret)[,Qval := p.adjust(Pval, method=padjust)]

  ret
}

#' @noRd
#' @import data.table
#' @importFrom dplyr filter group_by mutate summarize
.t.statisic.plate <- function(obj, mu)
{
  ret <- obj %>% ungroup
  mu <- .set.mu(ret, mu)
  ret <- dplyr::group_by(ret, Condition, Screen, Library,
                         ReadoutType, ScreenType, Cell, Design,
                         GeneSymbol, Entrez, Plate, Control,
                         RowIdx, ColIdx, Perturbation) %>%
    dplyr::summarize(Pval=.t.test(GeneSymbol, Readout, mu)$p.value,
                     Readout=mean(Readout, na.rm=TRUE))  %>%
    ungroup
  ret
}

#' @noRd
#' @import data.table
#' @importFrom dplyr filter
.set.mu <- function(ret, mu.gene)
{
  if (mu.gene != "0")
  {
    cont.tab <- dplyr::filter(ret, Control  == -1)
    if (tolower(mu.gene) == "scrambled")
      cont.tab <- dplyr::filter(cont.tab, tolower(GeneSymbol) == "scrambled")
    if (nrow(cont.tab) == 0)
      stop("No controls found for criteria!")
    if (nrow(cont.tab) < 3)
      stop("Less than three controls found better take mu=0!")
  mu <- cont.tab$Readout
  }
  else if (mu.gene == "0")
  {
    mu <- 0
  }
  else message(paste("\t..taking mu=", mu, ".",sep=""))
  if (!is.numeric(mu)) stop("Please provide a numeric mu!")

  mu
}

#' @noRd
#' @importFrom stats t.test
.t.test <- function(g, val, mu)
{
  tst <- list(p.value=1)
  if (length(val) < 3)
  {
    warning(paste("<3 values provided for" , g, " -> returning 1"))
  }
  else
  {
    tryCatch ({
      if (mu == 0)
      {
        tst <- stats::t.test(val, mu=0, alternative="two.sided")
      }
      else
      {
        tst <- stats::t.test(val, y=mu, alternative="two.sided")
      }
    }, warning = function(war)
    { warning(paste(war, " -> setting one for", g)); },
    error = function(err)
    { warning(paste(err, " -> setting one for", g)); })
  }

  tst
}

#' @noRd
#' @import data.table
#' @importFrom dplyr group_by summarize ungroup filter select mutate
#' @importFrom metap sumlog
.prioritize.tstatistic <- function(obj,
                                   hit.ratio=0.5,
                                   effect.size=0,
                                   pval.threshold=0.05,
                                   qval.threshold=1)
{
  res <- dplyr::group_by(obj, Condition, Screen, Library,
                         ScreenType, ReadoutType,
                         Design, Cell,
                         GeneSymbol, Entrez,
                         Plate, RowIdx, ColIdx, Perturbation) %>%
    dplyr::mutate(Hit=(Pval <= pval.threshold &
                  abs(Readout) >= effect.size)) %>%
    ungroup %>%
    dplyr::group_by(Condition, Screen, Library,
                    ReadoutType, ScreenType,
                    Design, Cell,
                    GeneSymbol, Entrez) %>%
    priorititize.statistic() %>%
    dplyr::filter(Pval <= pval.threshold & Qval <= qval.threshold)

  res
}
