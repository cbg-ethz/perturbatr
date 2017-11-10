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

#' Create model data for an LMM
#'
#' @export
#'
#' @import data.table
#'
#' @param obj  an object for which LMM model.data is created
#' @param drop  decide if genes that are not found in every virus should
#'  be dropped
#' @param ignore  ignore siRNAS that are only found \code{ignore} many times
#' @param weights  weights to set for the siRNAs
#' @param rel.mat.path  target-relation matrix (TODO)
#'
#' @return  returns an \code{knockdown.lmm.data} object
#'
#' @examples
#'  data(rnaiscreen)
#'  rnaiscreen <- preprocess(rnaiscreen, normalize="log")
#'
#'  lmm.model.data <- set.lmm.model.data(rnaiscreen)
set.lmm.model.data <- function(obj,
                               drop=TRUE,
                               ignore=1,
                               weights=NULL,
                               rel.mat.path=NULL)
{
  UseMethod("set.lmm.model.data")
}

#' @export
#' @method set.lmm.model.data knockdown.normalized.data
set.lmm.model.data.knockdown.normalized.data <- function(obj,
                                                        drop=TRUE,
                                                        ignore=1,
                                                        weights=NULL,
                                                        rel.mat.path=NULL)
{
  .set.lmm.matrix(obj, drop, ignore, weights, rel.mat.path)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
.set.lmm.matrix <- function(obj, drop, ignore, weights=NULL, rel.mat.path=NULL)
{
  lmm.mat <-
    dplyr::select(obj@.data, Entrez, GeneSymbol,
                  Virus, Readout, Control, Library,
                  Cell, ScreenType, Design, ReadoutType) %>%
    dplyr::filter(!is.na(GeneSymbol)) %>%
    # dont take positive controls since these are different between
    # the pathonges negative control on the other hand should be fine
    dplyr::filter(Control != 1)
  data.table::setDT(lmm.mat)[, Weight :=
                                .weights(lmm.mat, weights, rel.mat.path)]
  data.table::setDT(lmm.mat)[, VG := paste(Virus, GeneSymbol, sep=":")]
  #  remove library: not needed any more due to setting of weights
  data.table::setDT(lmm.mat)[, Library := NULL]

  lmm.mat$Entrez        <- as.integer(lmm.mat$Entrez)
  lmm.mat$Virus         <- as.factor(lmm.mat$Virus)
  lmm.mat$VG            <- as.factor(lmm.mat$VG)
  lmm.mat$Cell          <- as.factor(lmm.mat$Cell)
  lmm.mat$ReadoutType   <- as.factor(lmm.mat$ReadoutType)
  lmm.mat$ScreenType    <- as.factor(lmm.mat$ScreenType)
  lmm.mat$Design        <- as.factor(lmm.mat$Design)
  lmm.mat$GeneSymbol    <- as.factor(lmm.mat$GeneSymbol)
  lmm.mat$Weight        <- as.double(lmm.mat$Weight)
  lmm.mat$Control       <- as.integer(lmm.mat$Control)
  lmm.mat <-
    dplyr::filter(lmm.mat, !is.na(Readout)) %>%
    dplyr::group_by(VG) %>%
    dplyr::mutate(cnt=n()) %>%
    ungroup %>%
    dplyr::filter(cnt >= ignore) %>%
    dplyr::select(-cnt)

  if (drop)
  {
    vir.cnt <- .leuniq(lmm.mat$Virus)
    lmm.mat <- dplyr::group_by(lmm.mat, GeneSymbol) %>%
      # count if the genes are in all viruses
      # and compare if it matches the virus count
      dplyr::mutate(drop=(length(unique(Virus)) != vir.cnt)) %>%
      ungroup %>%
      dplyr::filter(!drop) %>%
      dplyr::select(-drop)
  }
  lmm.mat <- droplevels(lmm.mat)

  new("knockdown.lmm.data", .data=data.table::as.data.table(lmm.mat))

}

#' @noRd
#' @importFrom assertthat assert_that
.weights <- function(obj, weights, rel.mat.path)
{
  if (is.null(weights)) return(1)
  else if (!is.list(weights)) stop("Please give a list argument")
  els <- names(weights)
  ret <- rep(1, nrow(obj))
  for (el in els)
  {
    idxs <- switch(
      el,
      "pooled" = which(obj$Design == "pooled"),
      "single" = which(obj$Design == "single"),
      stop("Please provide 'single'/'pooled' list names for setting weights"))
    ret[idxs] <- as.numeric(weights[[el]])
    message(paste("Setting", length(idxs), el,
                  "well weights to:", as.numeric(weights[[el]])))
  }
  assertthat::assert_that(length(ret) == nrow(obj))
  ret
}
