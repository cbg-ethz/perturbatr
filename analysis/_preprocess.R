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


library(dplyr)
library(tibble)


source("_preprocess_background.R")
source("_preprocess_bscore.R")
source("_preprocess_loess.R")
source("_preprocess_log.R")
source("_preprocess_low_viability.R")
source("_preprocess_raw_data.R")
source("_preprocess_well_outliers.R")
source("_preprocess_zscore.R")


#' @title Pre-processing for perturbation data
#'
#' @description \code{preprocess} can be used for various preprocessing tasks,
#'  such as normalization, summarization, cell removal, removal of cytotoxic
#'  siRNAs or gRNAs.
#'
#'
#' @param obj  the object to be normalized
#' @param normalize  a vector with normalization methods
#' \itemize{
#'  \item{\emph{log} }{ does a log transformation}
#'  \item{\emph{poc} }{ divides by the mean/median of the control}
#'  \item{\emph{loess} }{ fits a local regression model to correct
#'                        for cell numbers per well}
#'  \item{\emph{b.score} }{ uses Tukey's median polish to correct
#'                          for row and column effects. B-scores can only be
#'                          used if data is provided in plate format, i.e.
#'                          when  \code{obj} contains columns named 'Plate',
#'                          'RowIdx' and 'ColIdx'. Otherwise an error will be
#'                           thrown.}
#'  \item{\emph{z.score} }{ standardizes a the data set.}
#'  \item{\emph{robust-z.score} }{ standardizes the data set using the median
#'                                 instead of mean and the MAD instead if the
#'                                 standard deviation,}
#'  \item{\emph{background} }{ substracts mean/median of some row/column or
#'                             all NA genes. This }
#'  \item{\emph{qq} }{ quantile-quantile normalization on replicates}
#' }
#' @param normalize.viability  boolean flag if the viability should also be
#'  normalized. This is only important if \code{obj} contains a column named
#'  'ReadoutClass', where some of the rows are labelled 'Readout' and some
#'  are labelled 'Viability'.
#' @param rm.cytotoxic  gene name of neutral gene against which viability is
#'  testet (e.g. \emph{Scrambled} or \emph{PIK3CA}), etc.
#' @param rm.outlier.wells  vector of quantiles, for instance `c(.05, .95)`
#'   The readout of wells with cell counts outside the quantile boundaries are
#'   set to NA.
#' @param z.score.level  if z-scores are used for normalisation, choose either
#'  'plate', if you want to correct on a plate-wise level and your data is
#'  in plate format, 'control' if you want to correct against a negative
#'
#' @param z.score.mu  the gene you want to compare to, e.g. 'scrambled', 'tp53'
#' @param poc.ctrl  the control used for poc-normalization, e.g. 'scrambled'
#' @param background.row  index of the row to correct backgorund against
#' @param background.column  index of column to correct backgorund against
#' @param summarization  method of how to summarize certain criteria
#'  (i.e. the background correction). Can be either 'mean' or 'median'.
#' @param drop  boolean if rows that are not needed should be dropped or
#'   kept (e.g. if you want to check whether all worked out correctly)
setGeneric(
  "preprocess",
  function(obj,
           normalize          = c("log", "robust-z.score"),
           normalize.viability= FALSE,
           rm.cytotoxic       = NULL,
           rm.outlier.wells   = NULL,
           z.score.level      = c("plate", "control"),
           z.score.mu         = "scrambled",
           poc.ctrl           = "scrambled",
           background.row     = NULL,
           background.column  = NULL,
           summarization      = c("mean", "median"),
           drop               = TRUE)
  {
    standardGeneric("preprocess")
  }
)


setMethod(
  "preprocess",
  signature = signature(obj="tbl_df"),
  function(obj,
          normalize           = c("log", "robust-z.score"),
          normalize.viability = FALSE,
          rm.cytotoxic        = NULL,
          rm.outlier.wells    = NULL,
          z.score.level       = c("plate", "control"),
          z.score.mu          = "scrambled",
          poc.ctrl            = "scrambled",
          background.row      = NULL,
          background.column   = NULL,
          summarization       = c("mean", "median"),
          drop                = TRUE)
  {
    res <- obj

    stopifnot(is.logical(drop))
    stopifnot(is.logical(normalize.viability))

    z.score.level    <- match.arg(z.score.level)
    summarization    <- match.arg(summarization)

    # do outlier removal based on cell numbers
    if (!is.null(rm.outlier.wells))
    {
      res <- .remove.outliers(res, rm.outlier.wells)
    }

    # do normalization
    res   <- .normalize(res,
                        normalize           = normalize,
                        normalize.viability = normalize.viability,
                        z.score.level       = z.score.level,
                        z.score.mu          = z.score.mu,
                        summarization       = summarization,
                        background.row      = background.row,
                        background.column   = background.column,
                        poc.ctrl            = poc.ctrl)

    # remove outliers based on cytotoxicity
    res <- .rm.cytotoxic(res, rm.cytotoxic)
    res <- .drop(res, drop)
    ret <- tibble::as.tibble(res)

    ret
  }
)


.drop <- function(obj, drop)
{
  if (drop)
  {
    if ("Viability" %in% colnames(obj))
      obj <- dplyr::select(obj, -Viability)
    obj <- dplyr::filter(obj, !is.na(GeneSymbol), GeneSymbol != "buffer")  %>%
      dplyr::select(-NumCells)
    if ("Remove" %in% colnames(obj))
      obj <- dplyr::filter(obj, Remove==FALSE) %>%
      dplyr::select(-Remove)
  }
  obj
}


.summarization.method <- function(summ.method)
{
  f <- switch(as.character(summ.method),
              "mean"=base::mean,
              "median"=stats::median,
              "min"=base::min,
              `NA`=NA,
              stop("wrong method given"))
  f
}
