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

#' Calculate statistics based on the chi-square-distribution to analyse the
#'  data.
#'
#' For this you should use a standardization before.
#'
#' @export
#' @import data.table
#'
#' @param obj  the data to be analysed
#' @param padjust  multiple testing correction method
#' @param ...   additional params
chisq.statistic <- function(obj, padjust=c("BH", "bonferroni"), ...)
{
  UseMethod("chisq.statistic", obj)
}

#' @export
#' @import data.table
#' @method chisq.statistic svd.data
chisq.statistic.svd.data <- function(obj, padjust=c("BH", "bonferroni"), ...)
{
  warning("this is a protoype. limited usability")
  stopifnot(length(maha) == nrow(obj))
  stopifnot(length(pvals) == nrow(obj))
  padjust <- match.arg(padjust)
  maha    <- .mahalanobis(obj$Readout)
  pvals   <- .chisq(maha)
  ret     <- obj
  data.table::setDT(ret)[, Pval := pvals]
  data.table::setDT(ret)[, Qval := p.adjust(pvals, method=padjust)]
  class(ret) <- c("svd.analysed.chisq", class(ret))
  ret
}

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
