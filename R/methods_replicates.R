# knockdown: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2015 - 2016 Simon Dirmeier
#
# This file is part of knockdown
#
# knockdown is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# knockdown is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with knockdown. If not, see <http://www.gnu.org/licenses/>.


#' Get the replicates of a data-set
#'
#' @export
#'
#' @param obj  an object for which the replicates should be retrieved
#' @param ...  additional params
#'
#' @return returns a \code{knockdown.replicates} object
#'
replicates <- function(obj, ...) UseMethod("replicates")

#' @export
#' @method replicates knockdown.raw.data
#' @import data.table
#' @importFrom dplyr filter
replicates.knockdown.raw.data <- function(obj, ...)
{
  obj@.data <- dplyr::filter(obj@.data, ReadoutClass=="Readout") %>%
    as.data.table
  replicates.knockdown.normalized.data(obj, ...)
}

#' @export
#' @method replicates knockdown.normalized.data
#' @import data.table
#' @importFrom dplyr filter
replicates.knockdown.normalized.data <- function(obj, ...)
{
  res <- obj@.data
  res <- dplyr::group_by(res, Virus, Screen, Replicate,
                            Cell, Design, Library,
                            ReadoutType, ScreenType) %>%
    dplyr::mutate(ReplicateIndex = .GRP) %>%
    ungroup

  new("knockdown.replicates", .data=as.data.table(res))
}
