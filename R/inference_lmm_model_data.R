# perturbatr: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2015 - 2016 Simon Dirmeier
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
# along with perturbatr If not, see <http://www.gnu.org/licenses/>.


#' Create model data for an hierarchical model
#'
#' @export
#'
#' @import data.table
#'
#' @param obj  an object for which hm model.data is created
#' @param drop  decide if genes that are not found in every Condition should
#'  be dropped
#' @param ignore  ignore siRNAS that are only found \code{ignore} many times
#' @param weights  weights to set for the siRNAs
#' @param rel.mat.path  target-relation matrix (TODO)
#'
#' @return  returns an \code{perturbation.hm.data} object
#' @examples
#'   data(rnaiscreen)
#'   rnaiscreen.normalized <- preprocess(rnaiscreen)
#'   set.hm.model.data(rnaiscreen.normalized)
set.hm.model.data <- function(obj,
                              drop=TRUE,
                              ignore=1,
                              weights=NULL,
                              rel.mat.path=NULL)
{
  UseMethod("set.hm.model.data")
}

#' @export
#' @method set.hm.model.data perturbation.normalized.data
set.hm.model.data.perturbation.normalized.data <- function(
  obj, drop=TRUE, ignore=1, weights=NULL, rel.mat.path=NULL)
{
  .set.hm.matrix(obj, drop, ignore, weights, rel.mat.path)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
.set.hm.matrix <- function(obj, drop, ignore, weights=NULL, rel.mat.path=NULL)
{
  hm.mat <-
    dplyr::select(obj@.data, Entrez, GeneSymbol,
                  Condition, Readout, Control, Library,
                  Cell, ScreenType, Design, ReadoutType) %>%
    dplyr::filter(!is.na(GeneSymbol)) %>%
    # dont take positive controls since these are different between
    # the pathonges negative control on the other hand should be fine
    dplyr::filter(Control != 1)
  data.table::setDT(hm.mat)[, Weight :=
                                .weights(hm.mat, weights, rel.mat.path)]
  data.table::setDT(hm.mat)[, VG := paste(Condition, GeneSymbol, sep=":")]
  #  remove library: not needed any more due to setting of weights
  data.table::setDT(hm.mat)[, Library := NULL]

  hm.mat$Entrez        <- as.integer(hm.mat$Entrez)
  hm.mat$Condition         <- as.factor(hm.mat$Condition)
  hm.mat$VG            <- as.factor(hm.mat$VG)
  hm.mat$Cell          <- as.factor(hm.mat$Cell)
  hm.mat$ReadoutType   <- as.factor(hm.mat$ReadoutType)
  hm.mat$ScreenType    <- as.factor(hm.mat$ScreenType)
  hm.mat$Design        <- as.factor(hm.mat$Design)
  hm.mat$GeneSymbol    <- as.factor(hm.mat$GeneSymbol)
  hm.mat$Weight        <- as.double(hm.mat$Weight)
  hm.mat$Control       <- as.integer(hm.mat$Control)

  hm.mat <-
    dplyr::filter(hm.mat, !is.na(Readout)) %>%
    dplyr::group_by(VG) %>%
    dplyr::mutate(cnt=n()) %>%
    ungroup %>%
    dplyr::filter(cnt >= ignore) %>%
    dplyr::select(-cnt)

  if (drop)
  {
    vir.cnt <- .leuniq(hm.mat$Condition)
    hm.mat <- dplyr::group_by(hm.mat, GeneSymbol) %>%
      dplyr::mutate(drop=(length(unique(Condition)) != vir.cnt)) %>%
      ungroup %>%
      dplyr::filter(!drop) %>%
      dplyr::select(-drop)
  }

  new("perturbation.hm.data", .data=data.table::as.data.table(droplevels(hm.mat)))

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
