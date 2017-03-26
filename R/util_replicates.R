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

#' Get the replicates of a data-set
#'
#' @export
#'
#' @param obj  an object for which the replicates should be retrieved
#' @param ...  additional params
replicates <- function(obj, ...) UseMethod("replicates")

#' @noRd
#' @export
#'
#' @import data.table
#' @importFrom dplyr filter
replicates.svd.raw <- function(obj, ...)
{
  obj <- dplyr::filter(obj, ReadoutClass=="Readout") %>% as.data.table
  replicates.svd.data(obj, ...)
}

#' @noRd
#' @export
#'
#' @import data.table
#' @importFrom dplyr filter
replicates.svd.data <- function(obj, ...)
{
  g <-
    dplyr::group_indices(obj, Virus, Screen, Replicate,
                         Design, Library,
                         ReadoutType, ScreenType) %>%
    as.data.table
  rep.frame <- obj
  rep.frame$grp <- g
  grps <- unique(rep.frame$grp)
  ret <- lapply(grps, function(i)
  {
    pl <- data.table::as.data.table(dplyr::filter(rep.frame, grp==i))
    class(pl) <- c("svd.replicate", class(pl))
    pl
  })
  class(ret) <- "svd.replicates"
  invisible(ret)
}
