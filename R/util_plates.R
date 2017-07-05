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


#' Get the plates of a data-set
#'
#' @export
#'
#' @param obj  an object for which the plates are going to be retrieved
#' @param ...  additional params
plates <- function(obj, ...)
{
  UseMethod("plates", obj)
}

#' @export
#' @method plates knockout.raw.data
#' @import data.table
#' @importFrom dplyr filter
plates.knockout.raw.data <- function(obj, ...)
{
  obj <- filter(obj, ReadoutClass=="Readout")
  plates.knockout.normalized.data(obj, ...)
}

#' @export
#' @method plates knockout.normalized.data
#' @import data.table
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
#' @importFrom dplyr group_indices
plates.knockout.normalized.data <- function(obj, ...)
{
  res <- obj@.data
  res   <- dplyr::group_by(res,
                         Virus, Screen, Replicate, Plate,
                         Cell, Design, Library,
                         ReadoutType, ScreenType) %>%
    dplyr::mutate(PlateIndex=.GRP) %>%
    ungroup

  new("knockout.plates", .data=as.data.table(res))
}
