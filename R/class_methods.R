# perturbR: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2018 Simon Dirmeier
#
# This file is part of perturbR
#
# perturbR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# perturbR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with perturbR. If not, see <http://www.gnu.org/licenses/>.


#' @include class_plates.R
#' @include class_replicates.R
#' @include class_quality.R


#' @title Extract parts of an object
#'
#' @description Extract a single plate from a \code{perturbation.plates} object
#'
#' @importFrom dplyr filter
#'
#' @param x  a \code{perturbation.plates} object
#' @param i  index of the plate to subset
#'
#' @return returns a \code{perturbation.plate} object
#'
setMethod(
  "[",
  signature=signature(x="perturbation.plates", i="numeric",
                      j="missing", drop="missing") ,
  function(x, i)
  {
    stopifnot(length(i) == 1)
    stopifnot(i >= 1)
    res <- x@.data
    if (max(res$PlateIndex) < i) stop("ArrayIndexOutOfBounds")
    res <- dplyr::filter(res, PlateIndex==i)

    new("perturbation.plate", .data=res)
  }
)

#' @title Extract parts of an object
#'
#' @description Extract a single replicate from a \code{perturbation.replicates} object
#'
#' @rdname subset-methods
#'
#' @importFrom dplyr filter
#'
#' @param x  a \code{perturbation.replicates} object
#' @param i  index of the plate to subset
#'
#' @return returns a \code{perturbation.replicate} object
#'
setMethod(
  "[",
  signature=signature(x="perturbation.replicates", i="numeric",
                      j="missing", drop="missing"),
  function(x, i)
  {
    stopifnot(length(i) == 1)
    stopifnot(i >= 1)
    res <- x@.data
    if (max(res$ReplicateIndex) < i) stop("ArrayIndexOutOfBounds")
    res <- dplyr::filter(res, ReplicateIndex==i)

    new("perturbation.replicate", .data=res)
  }
)
