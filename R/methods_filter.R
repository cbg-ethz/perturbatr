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


#' @title Filter the rows of a perturbation data set
#'
#' @description Takes a perturbation data set and filters the rows by some
#'  criterion. The filtered object will have the same type as the previous
#'  object.
#'
#' @export
#' @rdname filter-methods
#' @import data.table
#'
#' @param obj  the object to be filtered
#' @param ...  variable number of logical predicates in terms of the column
#'  names in \code{obj}. Multiple predicates will be combined with a logical
#'  'and'. Rows where all conditions are met will be kept. The column names do
#'   not need to be quoted. Wraps around
#'  \code{dplyr::filter}.
#'
#' @return  returns an object of the same type filtered by some criteria
#'
#' @examples
#'  data(rnaiscreen)
#'  flt.dat <- filter(rnaiscreen, Condition=="V1")
#'  flt.dat <- filter(rnaiscreen, Condition=="V1", RowIdx==1)
#'
filter <- function(obj, ...) UseMethod("filter")


#' @export
#' @method filter PerturbationData
#' @import data.table
#' @importFrom dplyr filter_
#' @importFrom lazyeval lazy_dots
#' @importFrom methods new
filter.PerturbationData <- function(obj, ...)
{
  filt.dat <- dplyr::filter_(dataSet(obj), .dots = lazyeval::lazy_dots(...))
  methods::new(class(obj)[1], dataSet=data.table::as.data.table(filt.dat))
}
