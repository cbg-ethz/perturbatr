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


#' @title Filter the rows of a knockout data set
#'
#' @description \code{lmm} TODO
#'
#' @export
#' @rdname filter-methods
#' @import data.table
#'
#' @param obj  the object to be filtered
#' @param ...  additional parameters
filter <- function(obj, ...) UseMethod("filter")


#' @export
#' @method filter knockout.data
#' @importFrom dplyr filter_
#' @importFrom lazyeval lazy_dots
filter.knockout.data <- function(obj, ...)
{
  filt.dat <- dplyr::filter_(obj@.data, .dots = lazyeval::lazy_dots(...))
  new(class(obj)[1], .data=filt.dat)
}
