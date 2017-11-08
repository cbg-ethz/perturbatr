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
#' @method bootstrap knockdown.data
#' @import data.table
#' @importFrom dplyr left_join mutate select group_by filter
#' @importFrom parallel mclapply detectCores
#' @importFrom methods new
bootstrap.knockdown.data <- function(obj, level=c("sirna", "pathogen"))
{
  dat <-
    dplyr::group_by(obj@.data, Virus, ScreenType, GeneSymbol) %>%
    dplyr::mutate(cnt=n(), grp=.GRP) %>%
    ungroup

  res  <- data.table::rbindlist(
    parallel::mclapply(
      unique(dat$grp),
      function (g)
      {
        grp.dat <- dplyr::filter(dat, grp==g)
        idx     <- sample(seq(grp.dat$cnt[1]), replace=TRUE)
        grp.dat[idx]
      },
      mc.cores=ifelse(tolower(Sys.info()['sysname']) %in% c("darwin", "unix"),
                      max(1, parallel::detectCores() - 1), 1L)
    )
  )

  ret <- methods::new(class(obj)[1], .data=dplyr::select(res, -cnt, -grp))
  ret
}

#' @noRd
loocv <- function(model.data, idx)
{
  UseMethod("loocv")
}

#' @method loocv svd.lmm.model.data
#' @import data.table
#' @importFrom dplyr left_join mutate select group_by filter
loocv.svd.lmm.model.data <- function(model.data, idx)
{
  if (!is.numeric(idx)) stop("provide an numeric index pls")
  dat <-
    model.data %>%
    dplyr::group_by(Virus, ScreenType, GeneSymbol) %>%
    dplyr::mutate(cnt=n(), grp=.GRP) %>%
    ungroup
  grps <- unique(dat$grp)
  res <- do.call(
    "rbind",
    lapply(
      grps,
      function (g)
      {
        grp.dat <- dplyr::filter(dat, grp==g)
        ma <- max(grp.dat$cnt)
        if (idx > ma || ma <= 3) {
          idxs <- seq(ma)
        } else {
          idxs <- seq(ma)[-idx]
        }
        grp.dat[idxs]
      }
    )
  )

  ret <- as.svd.lmm.model.data(res)
  ret
}
