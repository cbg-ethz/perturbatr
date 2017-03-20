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

#' Create model data for an LMM
#'
#' @export
#'
#' @import data.table
#'
#' @param obj  an object for which LMM model.data is created
#' @param drop  decide if genes that are not found in every virus should be dropped
#' @param ignore  ignore siRNAS that are only found \code{ignore} many times
#' @param weights  weights to set for the siRNAs
#' @param rel.mat.path  target-relation matrix (TODO)
set.lmm.model.data <- function(obj, drop=T, ignore=1, weights=NULL, rel.mat.path=NULL)
{
  UseMethod("set.lmm.model.data")
}

#' @export
#' @method set.lmm.model.data svd.data
set.lmm.model.data.svd.data <- function(obj, drop=T, ignore=1,
                                    weights=NULL, rel.mat.path=NULL)
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
  # setup pmm data
  pmm.mat <-
    # subset the columns
    dplyr::select(obj, Entrez, GeneSymbol, Virus, Readout, Control, Library,
                  Cell, ScreenType, Design, ReadoutType) %>%
    # only take entries that have a genesymbol
    dplyr::filter(!is.na(GeneSymbol)) %>%
    # dont take positive controls since these are different between the pathonges
    # negative control on the other hand should be fine
    dplyr::filter(Control != 1)
  # set weights
  data.table::setDT(pmm.mat)[, Weight := .weights(pmm.mat, weights, rel.mat.path)]
  # set a column that concats virus and genesymbol
  data.table::setDT(pmm.mat)[, VG := paste(Virus, GeneSymbol, sep=":")]
  #  remove librarz
  data.table::setDT(pmm.mat)[, Library := NULL]
  ## cast some columns
  pmm.mat$Entrez        <- as.integer(pmm.mat$Entrez)
  pmm.mat$Virus         <- as.factor(pmm.mat$Virus)
  pmm.mat$VG            <- as.factor(pmm.mat$VG)
  pmm.mat$Cell          <- as.factor(pmm.mat$Cell)
  pmm.mat$ReadoutType   <- as.factor(pmm.mat$ReadoutType)
  pmm.mat$ScreenType    <- as.factor(pmm.mat$ScreenType)
  pmm.mat$Design        <- as.factor(pmm.mat$Design)
  pmm.mat$GeneSymbol    <- as.factor(pmm.mat$GeneSymbol)
  pmm.mat$Weight        <- as.double(pmm.mat$Weight)
  pmm.mat$Control       <- as.integer(pmm.mat$Control)
  pmm.mat <-
    # remove entries that have nan as readout
    dplyr::filter(pmm.mat, !is.na(Readout)) %>%
    dplyr::group_by(VG) %>%
    # count how often VG is in the data
    dplyr::mutate(cnt=n()) %>%
    ungroup %>%
    # throw away VGs that are less than ignore
    dplyr::filter(cnt >= ignore) %>%
    # remove cnt column
    dplyr::select(-cnt)
  # drop genes that are not found in ALL pathogens
  if (drop)
  {
    # count how many viruses are in the date sets
    vir.cnt <- length(unique(pmm.mat$Virus))
    pmm.mat <-
      # group by genesymbol
      dplyr::group_by(pmm.mat, GeneSymbol) %>%
      # count if the genes are in all viruses
      # and compare if it matches the virus count
      dplyr::mutate(drop=(length(unique(Virus)) != vir.cnt)) %>%
      ungroup %>%
      # remove genes that are not in all genes
      dplyr::filter(!drop) %>%
      # remove drop column
      dplyr::select(-drop)
  }
  pmm.mat <- droplevels(pmm.mat)
  class(pmm.mat) <- c("svd.lmm.model.data", class(pmm.mat))
  pmm.mat
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
    idxs <- switch(el,
                   "pooled"  =which(obj$Design == "pooled"),
                   "single"=which(obj$Design == "single"),
                   stop("Please provide 'single'/'pooled' list names for setting weights"))
    ret[idxs] <- as.numeric(weights[[el]])
    message(paste("Setting", length(idxs), el,
                  "well weights to:", as.numeric(weights[[el]])))
  }
  assertthat::assert_that(length(ret) == nrow(obj))
  ret
}
