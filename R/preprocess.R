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


#' @include class-knockout_data.R


#' @title Pre-processing routine for data normalization, summarization and removal
#' of low viability siRNAs.
#'
#' @description \code{preproces} normalizes the data for comparable phenotypes.
#' TODO
#'
#' @export
#' @docType methods
#' @rdname preprocess-methods
#'
#' @import data.table
#'
#' @param obj  the object to be normalized
#' @param normalize  a vector with normalization methods
#' \itemize{
#'  \item{\emph{log} }{ does a log transformation}
#'  \item{\emph{poc} }{ divides by the mean/median of the control}
#'  \item{\emph{loess} }{ fits a local regression model to correct
#'                        for cell numbers per well}
#'  \item{\emph{b.score} }{ uses Tukey's median polish to correct
#'                          for row and column effects}
#'  \item{\emph{z.score} }{ standardizes a plate}
#'  \item{\emph{robust-z.score} }{ standardizes a plate using the median
#'                                 instead of mean and the MAD instead if the
#'                                 standard deviation}
#'  \item{\emph{background} }{ substracts mean/median of some row/column or
#'                             all NA genes}
#'  \item{\emph{qq} }{ quantile-quantile normalization on replicates}
#' }
#' @param normalize.viability  boolean flag if the viability should also be
#'  normalized
#' @param rm.cytotoxic  gene name of neutral gene against which viability is
#'  testet (e.g. \emph{Scrambled} or \emph{PIK3CA}), etc.
#' @param rm.outlier.wells  remove wells that have a extreme number of cells
#'  (outliers), e.g. by taking quantiles or absolute numbers
#' @param z.score.level  if z-scores are used for normalisation, choose either
#'   'plate' or 'control' to calculate statistics for
#' @param z.score.mu  the gene you want to compare to, e.g. 'scrambled', 'tp53'
#' @param poc.ctrl  the control used for poc-normalization, e.g. 'scrambled'
#' @param background.row  index of the row to correct backgorund against
#' @param background.column  index of column to correct backgorund against
#' @param outlier.well.range outlier cells will be removed based on this range
#' @param summarization  method of how to summarize screens (needed where? :/)
#' @param drop  boolean if rows that are not needed should be dropped or
#'   kept (e.g. if you want to check whether all worked out correctly)
setGeneric(
  "preprocess",
  function(obj,
           normalize          = c("log", "robust-z.score"),
           normalize.viability= F,
           rm.cytotoxic       = NULL,
           rm.outlier.wells   = c("none", "quantile", "absolute"),
           z.score.level      = c("plate", "control"),
           z.score.mu         = "scrambled",
           poc.ctrl           = "scrambled",
           background.row     = NULL,
           background.column  = NULL,
           outlier.well.range = c(.05, .95),
           summarization      = c("mean", "median"),
           drop               = T)
  {
    standardGeneric("preprocess")
  },
  package="knockout"
)


#' @rdname preprocess-methods
#' @aliases preprocess,knockout.data.raw-method
#'
#' @import data.table
#' @importFrom tidyr spread
#' @importFrom methods hasArg new
#' @importFrom dplyr group_by mutate filter select
setMethod(
  "preprocess",
  signature = signature(obj="knockout.data"),
  function(obj,
          normalize          = c("log", "robust-z.score"),
          normalize.viability= F,
          rm.cytotoxic       = NULL,
          rm.outlier.wells   = c("none", "quantile", "absolute"),
          z.score.level      = c("plate", "control"),
          z.score.mu         = "scrambled",
          poc.ctrl           = "scrambled",
          background.row     = NULL,
          background.column  = NULL,
          outlier.well.range = c(.05, .95),
          summarization      = c("mean", "median"),
          drop               = T)
  {
    if(obj@.type != .data.types()$RAW) stop("Data is already normalized.")

    res <- obj@.data

    stopifnot(is.logical(drop))
    stopifnot(is.logical(normalize.viability))

    rm.outlier.wells <- match.arg(rm.outlier.wells)
    z.score.level    <- match.arg(z.score.level)
    summarization    <- match.arg(summarization)

    # do outlier removal based on cell numbers
    if (rm.outlier.wells != "none")
      res <- .remove.outliers(res, rm.outlier.wells, outlier.well.range)
    # do normalization
    res   <- .normalize(res, normalize, normalize.viability)

    # remove outliers based on cytotoxicity
    res <- .rm.cytotoxic(res, rm.cytotoxic)
    res <- .drop(res, drop)
    res <- data.table::as.data.table(res)

    new("knockout.data", .data=res, .type="normalized")
  }
)

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
      obj <- dplyr::filter(obj, Remove==F) %>% dplyr::select(-Remove)
  }
  obj
}
