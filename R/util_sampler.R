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


#' Create a bootstrap sample from a data-set
#'
#' @export
#'
#' @param obj  the object which data should be bootstrapped
#' @param level  boostrap either on sirnas or on pathogens
#' @return returns an object with boostrapped data
#'
#' @examples
#' \dontrun{
#'   data(rnaiscreen)
#'   bootstrap(rnaiscreen)
#' }
bootstrap <- function(obj, level=c("sirna", "pathogen"))
{
  UseMethod("bootstrap")
}

#' @export
#' @method bootstrap data.table
#' @import data.table
#' @importFrom tibble as.tibble
#' @importFrom dplyr mutate select group_by filter group_indices
bootstrap.data.table <- function(obj, level=c("sirna", "pathogen"))
{
  dat <- tibble::as.tibble(obj)
  grps <- dplyr::group_indices(dat, Virus, ScreenType, GeneSymbol)
  dat  <- dplyr::mutate(dat, grp=grps) %>%
    dplyr::group_by(Virus, ScreenType, GeneSymbol) %>%
    dplyr::mutate(cnt=n()) %>%
    ungroup

  res <- do.call(
      "rbind",
      lapply(unique(dat$grp),
      function (g)
      {
        grp.dat <- dplyr::filter(dat, grp==g)
        idx     <- sample(seq(grp.dat$cnt[1]), replace=TRUE)
        grp.dat[idx,]
      }))

  ret <- dplyr::select(res, -cnt, -grp)
  ret
}


#' @export
#' @method bootstrap knockdown.data
#' @import data.table
#' @importFrom methods new
bootstrap.knockdown.data <- function(obj, level=c("sirna", "pathogen"))
{
  res <- bootstrap(obj@.data)
  ret <- methods::new(class(obj)[1], .data=data.table::as.data.table(res))
  ret
}
