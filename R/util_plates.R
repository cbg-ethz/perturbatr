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
  g <- dplyr::group_indices(obj@.data, Virus, Screen, Replicate, Plate,
                            ReadoutType, ScreenType) %>%
    as.data.table
  plate.frame     <- obj
  plate.frame$grp <- g
  grps <- unique(plate.frame$grp)
  plates <- lapply(grps, function(i)
  {
    pl <- data.table::as.data.table(dplyr::filter(plate.frame, grp==i))
    class(pl) <- c("svd.plate", class(pl))
    pl
  })
  class(plates) <- "svd.plates"
  invisible(plates)
}

#' @noRd
#' @export
#' @method plates default
plates.default <- function(obj, ...)
{
  # TODO make better
  plate.frame <- obj
  class(plate.frame) <- c("svd.plate.rows", "svd.plate", class(plate.frame))
  plate.frame
}
