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


#' @title Calculate statistics based on the chi-square-distribution to analyse the
#'  data.
#'
#' @description TODO
#'
#' @export
#' @docType methods
#' @rdname chisq_statistic-methods
#'
#' @param obj  the data to be analysed
#' @param padjust  multiple testing correction method
#' @param ...   additional params
setGeneric(
  "chisq.statistic",
  function(obj, padjust=c("BH", "bonferroni"), ...)
  {
      standardGeneric("chisq.statistic")
  },
  package="knockout"
)

#' @rdname chisq_statistic-methods
#' @aliases chisq.statistic,knockout.lmm.data-method
#' @import data.table
setMethod(
  "chisq.statistic",
  signature = signature(obj="knockout.data"),
  function(obj, padjust=c("BH", "bonferroni"), ...)
  {
    .check.data(obj)
    dat <- obj@.data
    if (.leuniq(dat$Replicate) !=  1)
      stop(paste0("You provided a data-set with several replicates. ",
                  "Summarize these first or use another method."))
    padjust <- match.arg(padjust)
    maha    <- .mahalanobis(dat$Readout)
    pvals   <- .chisq(maha)
    stopifnot(length(maha) == nrow(dat), length(pvals) == nrow(dat))
    data.table::setDT(dat)[, Effect := maha]
    data.table::setDT(dat)[, Pval   := pvals]
    data.table::setDT(dat)[, Qval   := p.adjust(pvals, method=padjust)]

    ret <- new("knockout.analysed",
               .inference=.inference.types()$CHISQ.TEST,
               .data=dat)

  }
)

#' @export
#' @method chisq.statistic numeric
#' @importFrom stats p.adjust
chisq.statistic.numeric <- function(obj, padjust=c("BH", "bonferroni"), ...)
{
  warning("this is a protoype. limited usability")
  padjust  <- match.arg(padjust)
  maha     <- .mahalanobis(obj)
  p.vals   <- .chisq(maha)
  list(p.vals=p.vals, q.vals=stats::p.adjust(p.vals, method=padjust))
}

#' @noRd
#' @importFrom stats sd
.mahalanobis <- function(vals)
{
  sd   <- stats::sd(vals, na.rm=T)
  mu   <- mean(vals, na.rm=T)
  rm   <- (vals - mu)
  maha <- sqrt( rm * rm / sd)
  maha
}

#' @noRd
#' @importFrom stats pchisq
.chisq <- function(vals) pchisq(vals, df=1, lower.tail=F)
