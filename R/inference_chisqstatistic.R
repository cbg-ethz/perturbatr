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


#' @include class_knockout_data.R


#' @title Calculate statistics based on the chi-square-distribution to analyse
#'  the data.
#'
#' @description TODO
#'
#' @export
#' @docType methods
#' @rdname chisq_statistic-methods
#'
#' @param obj  the data to be analysed
#' @param padjust  multiple testing correction method
#' @param effect.size  the relative strength of a signal to count as a hit,
#'  i. e. the biological significance of a gene/sirna
#' @param pval.threshold  the significance level for a sirna/gene,
#'  i. e. the statistical significance of a gene/sirna
#'  to be counted as significant
#' @param qval.threshold  the significance level of the multiple testing
#'  corrected p-value. This should be set to an appropriate significance level
#'  just like the \code{pval.threshold} as well.
#'
#' @return returns a \code{knockout.chisqstatistic.analysed} object
#'
#' @examples
#'  data(rnaiscreen)
#'  screen.norm <- preprocess(rnaiscreen, normalize="log")
#'  screen.norm <- filter(screen.norm, Replicate == 1)
#'
#'  res <- chisq.statistic(screen.norm)
setGeneric(
  "chisq.statistic",
  function(obj,
           padjust=c("BH", "bonferroni"),
           effect.size=0,
           pval.threshold=0.05,
           qval.threshold=1)
  {
      standardGeneric("chisq.statistic")
  },
  package="knockout"
)

#' @rdname chisq_statistic-methods
#' @aliases chisq.statistic,knockout.data-method
#' @import data.table
setMethod(
  "chisq.statistic",
  signature = signature(obj="knockout.data"),
  function(obj,
           padjust=c("BH", "bonferroni"),
           effect.size=0,
           pval.threshold=0.05,
           qval.threshold=1)
  {
    if (.leuniq(obj@.data$Replicate) !=  1)
      stop(paste0("You provided a data-set with several replicates. ",
                  "Summarize these first or use another method."))

    res <- .chisq.statistic(obj@.data, match.arg(padjust))
    priorit <- .prioritize.chisq.tstatistic(
      res, effect.size, pval.threshold, qval.threshold)

    ret <- methods::new(
      "knockout.chisqstatistic.analysed",
      .data = res,
      .gene.hits = data.table::as.data.table(priorit),
      .params = list(effect.size=effect.size,
                     pval.threshold=pval.threshold,
                     qval.threshold=qval.threshold)
    )

    ret
  }
)

#' @noRd
#' @import data.table
#' @importFrom dplyr mutate group_indices filter group_by
.chisq.statistic <- function(obj, padjust)
{
  ret <- dplyr::group_by(obj, Virus, Screen, ReadoutType,
                         ScreenType, Library, Design, Cell) %>%
    dplyr::mutate(grp = .GRP)
  grps <- unique(ret$grp)

  ret  <- do.call(
    "rbind",
    lapply
    (
      grps,
      function(g)
      {
        grp.dat <- dplyr::filter(ret, grp == g)
        message(paste("Doing grp: ", paste(grp.dat$Virus[1],
                                           grp.dat$Screen[1],
                                           grp.dat$ScreenType[1],
                                           grp.dat$ReadoutType[1],
                                           grp.dat$Cell[1],
                                           grp.dat$Design[1],
                                           grp.dat$Library[1],
                                           sep=", ")))
        fr <- .do.chisq(grp.dat, padjust)
        fr %>% dplyr::select(-grp)
      }
    )
  )

  ret
}

#' @noRd
#' @import data.table
#' @importFrom assertthat assert_that
.do.chisq <- function(obj, padjust)
{
  maha    <- .mahalanobis(obj$Readout)
  pvals   <- .chisq(maha)
  assertthat::assert_that(length(maha) == nrow(obj),
                          length(pvals) == nrow(obj))
  data.table::setDT(obj)[, Effect := maha]
  data.table::setDT(obj)[, Pval   := pvals]
  data.table::setDT(obj)[, Qval   := p.adjust(Pval, method=padjust)]

  obj
}

#' @noRd
#' @importFrom stats sd
.mahalanobis <- function(vals)
{
  sd   <- stats::sd(vals, na.rm=TRUE)
  mu   <- mean(vals, na.rm=TRUE)
  rm   <- (vals - mu)
  maha <- sqrt(rm * rm / sd)
  maha
}

#' @noRd
#' @importFrom stats pchisq
.chisq <- function(vals) pchisq(vals, df=1, lower.tail=FALSE)


#' @noRd
#' @import data.table
#' @importFrom dplyr group_by summarize ungroup filter select mutate
.prioritize.chisq.tstatistic <- function(obj,
                                         effect.size,
                                         pval.threshold,
                                         qval.threshold)
{
  res <- dplyr::filter(obj,
                       abs(Readout) >= effect.size,
                       Pval <= pval.threshold,
                       Qval <= qval.threshold)

  res
}
