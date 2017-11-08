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


#' Get the plates of a data-set
#'
#' @export
#'
#' @param obj  an object for which the plates are going to be retrieved
#' @param ...  additional params
#'
#' @return returns a \code{knockdown.plates} object
#'
#' @examples
#'  data(rnaiscreen)
#'  plates <- plates(rnaiscreen)
plates <- function(obj, ...)
{
  UseMethod("plates", obj)
}

#' @export
#' @method plates knockdown.raw.data
#' @import data.table
#' @importFrom dplyr filter
plates.knockdown.raw.data <- function(obj, ...)
{
  obj <- filter(obj, ReadoutClass=="Readout")
  plates.knockdown.normalized.data(obj, ...)
}

#' @export
#' @method plates knockdown.normalized.data
#' @import data.table
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
#' @importFrom dplyr group_indices
plates.knockdown.normalized.data <- function(obj, ...)
{
  res <- obj@.data
  res   <- dplyr::group_by(res,
                         Virus, Screen, Replicate, Plate,
                         Cell, Design, Library,
                         ReadoutType, ScreenType) %>%
    dplyr::mutate(PlateIndex=.GRP) %>%
    ungroup

  new("knockdown.plates", .data=as.data.table(res))
}
