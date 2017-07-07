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


#' @title Select columns of a knockout data set
#'
#' @description Takes a knockout data set and selects columns by some
#'  criterion. The filtered object will have class \code{data.table}.
#'
#' @export
#' @import data.table
#'
#' @param obj  the object of which columns should be selected
#' @param ...  additional parameters
#'
#' @return  returns a \code{data.table} with specified column
#'
#' @examples
#'  data(rnaiscreen)
#'  qual <- select(rnaiscreen, Virus)
select <- function(obj, ...) UseMethod("select")


#' @export
#' @method select knockout.data
#' @import data.table
#' @importFrom dplyr select_
#' @importFrom lazyeval lazy_dots
select.knockout.data <- function(obj, ...)
{
  dplyr::select_(obj@.data, .dots = lazyeval::lazy_dots(...)) %>%
    as.data.table
}
