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


#' @title Select columns of a perturbation data set
#'
#' @description Takes a perturbation data set and selects columns by some
#'  criterion. The filtered object will have class \code{data.table}.
#'
#' @export
#' @import data.table
#'
#' @param obj  the object of which columns should be selected
#' @param ...  variable number of unquoted columns names from \code{obj} to
#'  subset from. Negative column names will be deselected. Wraps around
#'  \code{dplyr::select}.
#'
#' @return  returns a \code{tibble} with the specified columns
#'
#' @examples
#'  data(rnaiscreen)
#'  perturbatr::select(rnaiscreen, GeneSymbol, Readout)
select <- function(obj, ...) UseMethod("select")


#' @export
#' @method select PerturbationData
#' @import data.table
#' @importFrom dplyr select_
#' @importFrom lazyeval lazy_dots
#' @importFrom tibble as_tibble
select.PerturbationData <- function(obj, ...)
{
  df <- dplyr::select_(dataSet(obj), .dots = lazyeval::lazy_dots(...))
  tibble::as.tibble(df)
}
