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


#' Extract parts of an object
#'
#' @description Extract a single plate from a \code{knockout.plates} object
#'
#' @importFrom dplyr filter
#'
#' @param x  a \code{knockout.plates} object
#' @param i  index of the plate to subset
#'
#' @return returns a \code{knockout.plate} object
#'
setMethod(
  "[",
  signature=signature(x="knockout.plates", i="numeric",
                      j="missing", drop="missing") ,
  function(x, i)
  {
    stopifnot(length(i) == 1)
    stopifnot(i >= 1)
    res <- x@.data
    if (max(res$PlateIndex) < i) stop("ArrayIndexOutOfBounds")
    res <- dplyr::filter(res, PlateIndex==i)

    new("knockout.plate", .data=res)
  }
)

#' Extract parts of an object
#'
#' Extract a single replicate from a \code{knockout.replicates} object
#'
#' @rdname subset-methods
#'
#' @importFrom dplyr filter
#'
#' @param x  a \code{knockout.replicates} object
#' @param i  index of the plate to subset
#'
#' @return returns a \code{knockout.replicate} object
#'
setMethod(
  "[",
  signature=signature(x="knockout.replicates", i="numeric",
                      j="missing", drop="missing"),
  function(x, i)
  {
    stopifnot(length(i) == 1)
    stopifnot(i >= 1)
    res <- x@.data
    if (max(res$ReplicateIndex) < i) stop("ArrayIndexOutOfBounds")
    res <- dplyr::filter(res, ReplicateIndex==i)

    new("knockout.replicate", .data=res)
  }
)
