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


#' @title Calculate statistics based on the t-test to analyse the data.
#'
#' @description For this you should use a standardization before.
#'
#' @export
#' @docType methods
#' @rdname t_statistic-methods
#'
#' @import data.table
#'
#' @param obj  the data to be analysed
#' @param mu  side to which the mean of a siRNA is compared to
#' @param padjust  multiple testing correction method
#' @param hit.ratio  the ratio of siRNAs
#' @param effect.size  the relative strength of a signal to count as a hit,
#'  i. e. the biological significance of a gene/sirna
#' @param pval.threshold  the significance level for a sirna/gene,
#'  i. e. the statistical significance of a gene/sirna
#'  to be counted as significant
#' @param qval.threshold  the significance level of the multiple testing
#'  corrected p-value. This should be set to an appropriate significance level
#'  just like the \code{pval.threshold} as well.
#'
#' @return returns a \code{knockdown.tstatistic.analysed} object
#'
#' @examples
#'  data(rnaiscreen)
#'  screen.norm <- preprocess(rnaiscreen, normalize="log")
#'
#'  res <- tstatistic(screen.norm)
#'
setGeneric(
  "tstatistic",
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
  package="knockdown"
)


#' @rdname t_statistic-methods
#' @aliases tstatistic,knockdown.data-method
#' @import data.table
#' @importFrom methods new
setMethod(
  "tstatistic",
  signature = signature(obj = "knockdown.data"),
  function(obj,
           mu=c(0, "scrambled", "control"),
           padjust=c("BH", "bonferroni"),
           hit.ratio=0.5,
           effect.size=0,
           pval.threshold=0.05,
           qval.threshold=1)
  {
    res     <- .t.statistic(obj@.data,
                            mu=match.arg(mu),
                            padjust=match.arg(padjust))
    priorit <- .prioritize.tstatistic(
      res, hit.ratio, effect.size, pval.threshold, qval.threshold)

    ret <- methods::new(
      "knockdown.tstatistic.analysed",
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
  ret <- dplyr::group_by(obj, Virus, Screen, Library,
                         ReadoutType, ScreenType,
                         Design, Cell) %>%
    dplyr::mutate(grp=.GRP) %>%
    ungroup
  grps <- unique(ret$grp)
  res <- do.call(
    "rbind",
    lapply(
      grps,
      function (g)
      {
        grp.dat <- dplyr::filter(ret, grp==g) %>% ungroup
        message(paste("Doing grp: ", paste(grp.dat$Virus[1],
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
      grps,
      function(g)
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
  ret <- dplyr::group_by(ret, Virus, Screen, Library,
                         ReadoutType, ScreenType, Cell, Design,
                         GeneSymbol, Entrez, Plate, Control,
                         RowIdx, ColIdx, siRNAIDs) %>%
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
  res <- dplyr::group_by(obj, Virus, Screen, Library,
                         ScreenType, ReadoutType,
                         Design, Cell,
                         GeneSymbol, Entrez,
                         Plate, RowIdx, ColIdx, siRNAIDs) %>%
    dplyr::mutate(Hit=(Pval <= pval.threshold &
                  abs(Readout) >= effect.size)) %>%
    ungroup %>%
    dplyr::group_by(Virus, Screen, Library,
                    ReadoutType, ScreenType,
                    Design, Cell,
                    GeneSymbol, Entrez) %>%
    dplyr::summarize(HitRatio   = (sum(Hit == TRUE, na.rm=TRUE)/n()),
                     Pval       = metap::sumlog(Pval)$p,
                     Qval       = metap::sumlog(Qval)$p,
                     MeanEffect = mean(Readout, na.rm=TRUE),
                     MaxEffect  = max(Readout, na.rm=TRUE),
                     MinEffect  = min(Readout, na.rm=TRUE),
                     AllPval=paste(sprintf("%03f", Pval), collapse=","),
                     AllQval=paste(sprintf("%03f", Qval), collapse=",")) %>%
    ungroup %>%
    dplyr::filter(Pval <= pval.threshold & Qval <= qval.threshold)

  res
}
