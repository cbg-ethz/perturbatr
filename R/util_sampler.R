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

#' @noRd
bootstrap <- function(model.data, level=c("sirna", "pathogen"))
{
  UseMethod("bootstrap")
}

#' @import data.table
#' @importFrom dplyr left_join mutate select group_by filter
#' @method bootstrap svd.lmm.model.data
bootstrap.svd.lmm.model.data <- function(model.data,
                                         level=c("sirna", "pathogen"))
{
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
        idx <- sample(seq(grp.dat$cnt[1]), replace=T) %>% unique
        grp.dat[idx]
      }
    )
  )
  res <- as.svd.lmm.model.data(res)
  res
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
  res <- as.svd.lmm.model.data(res)
  res
}
