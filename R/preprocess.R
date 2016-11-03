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

#' Pre-processing routine for data normalization, summarization and removal of low viability siRNAs.
#'
#' @export
#'
#' @import data.table
#'
#' @param obj  the data.table to be analysed
#' @param normalize  a vector with normalization methods
#' \itemize{
#'  \item{\emph{log} }{ does a log transformation}
#'  \item{\emph{poc} }{ divides by the mean/median of the control}
#'  \item{\emph{loess} }{ fits a local regression model to correct for cell numbers per well}
#'  \item{\emph{b.score} }{ uses Tukey's median polish to correct for row and column effects}
#'  \item{\emph{z.score} }{ standardizes a plate}
#'  \item{\emph{robust-z.score} }{ standardizes a plate using the median instead of mean and the MAD instead if the standard deviation}
#'  \item{\emph{background} }{ substracts mean/median of some row/column or all NA genes}
#'  \item{\emph{qq} }{ quantile-quantile normalization on replicates}
#' }
#' @param normalize.viability  boolean flag if the viability should also be removed
#' @param rm.cytotoxic  gene name of neutral gene against which viability is testet (e.g. \emph{Scrambled} or \emph{PIK3CA})
#' @param rm.outlier.wells  remove wells that have a extreme number of cells (outliers)
#' @param drop  boolean if rows that are not needed should be dropped
#' @param ...  additional parameters
#' \itemize{
#'  \item{\emph{method} }{ either \emph{mean} or \emph{median} for diverse summarizations/..}
#'  \item{\emph{drop} }{ boolean flag if useless columns should be dropped or not}
#'  \item{\emph{z.score.level} }{ either \emph{plate} or  \emph{control}, depending on which the z-scores should be calculated}
#'  \item{\emph{z.score.ctrl} }{ gene name for which plate is normalized to when level=control}
#'  \item{\emph{poc.ctrl} }{ gene name for which plate is normalized to when poc is used}
#'  \item{\emph{background.row} }{ row idx that is used for background correction}
#'  \item{\emph{background.column} }{ column idx that is used for background correction}
#'  \item{\emph{outlier.well.range} }{ vector of lower and upper bound for excluded cells (e.g.: c(150,500) or c(.05, .95) depending on your choice of rm.outlier.wells) }
#' }
preprocess <- function(obj, normalize=c("log", "robust-z.score"),
                       normalize.viability=F, rm.cytotoxic=NULL,
                       rm.outlier.wells=c(NA, "quantile"), drop=T, ...)
{
  UseMethod("preprocess", obj)
}

#' @noRd
#' @export
#' @import data.table
#' @importFrom dplyr group_by mutate filter select
#' @importFrom tidyr spread
preprocess.svd.raw <- function (obj,normalize=c("log", "robust-z.score"),
                                normalize.viability=F, rm.cytotoxic=NULL,
                                rm.outlier.wells=c(NA, "quantile"), drop=T, ...)
{
  params <- list(...)
  # what wells should set to NA dependent on quantiles
  outlier.well.range <- NA
  if (hasArg(outlier.well.range))
    outlier.well.range <- params$outlier.well.range
  rm.outlier.wells <- match.arg(rm.outlier.wells)
  if (!is.logical(drop))
    stop("Please provide a boolean for 'drop'!")
  # should viabilityies also be normalized
  if (!is.logical(normalize.viability))
    stop("Please provide a boolean normalize.viability")
  # do outlier removal
  res <- .remove.outliers(obj, rm.outlier.wells, outlier.well.range)
  # do normalization
  res <- .normalize(res, normalize, normalize.viability, ...)
  # remove outliers
  res <- .rm.cytotoxic(res, rm.cytotoxic)
  res <- .drop(res, drop)
  res <- as.data.table(res)
  class(res) <- c("svd.data", class(res))
  invisible(res)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr select filter
.drop <- function(obj, drop)
{
  if (drop)
  {
    if ("Viability" %in% colnames(obj))
      obj <- dplyr::select(obj, -Viability)
    obj <- dplyr::filter(obj, !is.na(GeneSymbol), GeneSymbol != "buffer")  %>%
      dplyr::select(-NumCells)
    if ("Remove" %in% colnames(obj))
      obj <- dplyr::select(obj, -Remove)
  }
  obj
}
