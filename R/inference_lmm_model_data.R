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
#' @docType methods
#' @rdname setModelData-methods
#'
#' @import data.table
#'
#' @param obj  an data set
#' @param drop  boolean if genes that are not found in every Condition should
#'  be dropped
#' @param weights  a numeric vector used as weights for the single
#'  perturbations
#'
#' @import data.table
#' @importFrom dplyr select filter group_by mutate
#'
#' @return  returns an \code{PerturbationData} object
#' @examples
#'   data(rnaiscreen)
#'   rnaiscreen.normalized <- preprocess(rnaiscreen)
#'   setModelData(rnaiscreen.normalized)
setGeneric(
  "setModelData",
  function(obj, drop=TRUE,  weights=NULL)
)


#' @rdname setModelData-methods
#' @aliases setModelData,PerturbationData-method
setMethod(
  "setModelData",
  signature = signature(obj="PerturbationData")
  function(obj, drop=TRUE, ignore=1, weights=NULL)
  {
    hm.mat <- dataSet(obj@) %>%
      dplyr::mutate(Weight = as.double(weights)) %>%
      dplyr::filter(!is.na(GeneSymbol)) %>%
      dplyr::filter(Control != 1) %>%
      dplyr::filter(!is.na(Readout))

    data.table::setDT(hm.mat)[, CG := paste(Condition, GeneSymbol, sep=":")]

    hm.mat$Condition     <- as.factor(hm.mat$Condition)
    hm.mat$CG            <- as.factor(hm.mat$CG)
    hm.mat$GeneSymbol    <- as.factor(hm.mat$GeneSymbol)

    if (drop)
    {
      vir.cnt <- .leuniq(hm.mat$Condition)
      hm.mat <- dplyr::group_by(hm.mat, GeneSymbol) %>%
        dplyr::mutate(drop=(length(unique(Condition)) != vir.cnt)) %>%
        ungroup %>%
        dplyr::filter(!drop) %>%
        dplyr::select(-drop)
    }

    hm.mat
  }
)
