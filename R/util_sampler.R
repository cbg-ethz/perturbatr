# perturbatr: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2018 Simon Dirmeier
#
# This file is part of perturbatr
#
# perturbatr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# perturbatr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with perturbatr. If not, see <http://www.gnu.org/licenses/>.


#' Create a bootstrap sample from a data-set
#'
#' @export
#'
#' @param obj  the object which data should be bootstrapped
#' @param ...  groups on which you be bootstrapped. If you want to
#'  create a normal boostrap sample, you would ignore this argument. If you
#'  want to separate your data into groups and bootstrap from every group, you
#'  would give the unquoted name of the columns in your \code{obj} to group on.
#'  For
#'  instance, if you provide `gene` as an argument, then your data set would be
#'  grouped into separate `gene` groups and bootstrapping would be conducted on
#'  every of those groups. Afterwards genes are aggregated
#'
#' @return returns an object with boostrapped data
#'
#' @examples
#'   data(rnaiscreen)
#'   bootstrap(rnaiscreen)
#'   bootstrap(rnaiscreen, Condition, Perturbation)
bootstrap <- function(obj, ...)
{
    UseMethod("bootstrap")
}

#' @export
#' @method bootstrap tbl_df
#' @import tibble
#' @import dplyr
bootstrap.tbl_df <- function(obj, ...)
{
  dat  <- tibble::as.tibble(obj)
  dat  <- dplyr::group_by_(dat, .dots = lazyeval::lazy_dots(...)) %>%
    { dplyr::mutate(dplyr::ungroup(.), grp = dplyr::group_indices(.)) } %>%
    dplyr::group_by_(.dots = lazyeval::lazy_dots(...)) %>%
    dplyr::mutate(cnt = n()) %>%
    dplyr::ungroup()

  res <- dplyr::bind_rows(
    lapply(unique(dat$grp),
    function (g)
    {
      grp.dat <- dplyr::filter(dat, grp == g)
      idx     <- sample(seq(grp.dat$cnt[1]), replace=TRUE)
      grp.dat[idx, ]
    }))

  ret <- dplyr::select(res, -cnt, -grp)
  ret
}


#' @export
#' @method bootstrap PerturbationData
#' @import tibble
#' @importFrom methods new
bootstrap.PerturbationData <- function(obj, ...)
{
  res <- bootstrap(dataSet(obj), ...)
  ret <- methods::new(class(obj)[1], dataSet = tibble::tibble(res))
  ret
}
