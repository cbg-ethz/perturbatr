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

#' Converts an object to an svd.data object
#'
#' @export
#' @import data.table
#' @param obj  the object to be converted
#' @param ...  additional params
as.svd.data <- function(obj, ...) UseMethod("as.svd.data")

#' @export
#' @method as.svd.data data.table
as.svd.data.data.table <- function(obj, ...)
{
  class(obj) <- unique(c("svd.data", class(obj)))
  invisible(obj)
}

#' Converts an object to an svd.lmm.model.data object
#'
#' @export
#' @import data.table
#' @param obj  the object to be converted
#' @param ...  additional params
as.svd.lmm.model.data <- function(obj, ...) UseMethod("as.svd.lmm.model.data")

#' @export
#' @method as.svd.lmm.model.data data.table
as.svd.lmm.model.data.data.table <- function(obj, ...)
{
  class(obj) <- unique(c("svd.lmm.model.data", class(obj)))
  invisible(obj)
}

#' Convert to an plate object
#'
#' @export
#'
#' @param obj  objject to be converted
#' @param ...  additional params
as.svd.plate <- function(obj, ...) UseMethod("as.svd.plate")
