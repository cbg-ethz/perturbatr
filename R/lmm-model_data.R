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
  UseMethod("model.data.lmm")
}

#' @export
#' @method set.lmm.model.data svd.data
set.lmm.model.data.svd.data <- function(obj, drop=T, ignore=1,
                                    weights=NULL, rel.mat.path=NULL)
{
  .set.lmm.matrix(obj, drop, ignore, weights, rel.mat.path)
}

#' @export
#' @method model.data.lmm data.table
model.data.lmm.data.table <- function(obj, drop=T, ignore=1,
                                      weights=NULL, rel.mat.path=NULL)
{
  .set.lmm.matrix(obj, drop, ignore, weights, rel.mat.path)
}
