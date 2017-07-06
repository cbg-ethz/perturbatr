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
#' @import data.table
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
.qq.norm <- function(obj, level=c("replicate", "plate"))
{
  stop("Not implemented yet")
  level <- match.arg(level)
  # TODO implement on replicates
  message(paste("Calculating quantile-quantile normalisation on", level))
  # TODO difference between plate/replicates
  ret <-
    dplyr::group_by(obj, Virus, Library, Screen,
                    ReadoutType, ScreenType, ReadoutClass,
                    Design, Cell) %>%
    dplyr::mutate(Readout = .qq.norm.group(Replicate, Readout,
                                           RowIdx, ColIdx)) %>%
    ungroup
  invisible(ret)
}

#' @noRd
#' @import data.table
#' @importFrom limma normalizeBetweenArrays
.qq.norm.group <- function(grps, readouts, row.idxs, col.idxs)
{
  if (length(unique(grps)) == 1)
    stop("Please provide multiple replicates/plates per screen!")
  grps.sizes    <- unname(table(grps))
  if (length(unique(grps.sizes)) != 1)
    stop("Please provide data with equal number of samples in every group.")
  ord <- data.table::data.table(Grp=paste("g", grps, sep=""),
                                Readouts=readouts,
                                Row=row.idxs, Col=col.idxs) %>%
    dplyr::group_by(Grp) %>%
    dplyr::mutate(ord=1:n()) %>%
    ungroup %>%
    tidyr::spread(Grp, Readouts)
  qq.mat <- dplyr::select(ord, -Row, -Col, -ord) %>% as.matrix
  res <- ord
  tryCatch ({
    norm.mat <- limma::normalizeBetweenArrays(qq.mat, method="quantile")
    res <- dplyr::select(ord, Row, Col, ord) %>% cbind(norm.mat) %>%
      as.data.table
  }, warning = function(war) { print(paste("WARNING: ", war)); },
  error = function(err) { print(paste("ERROR: ", err)); }
  )
  ret <- tidyr::gather(res, Grp, Readout, g1:g9) %>%
    as.data.table %>%
    .[, .SD[order(ord)], by="Grp"]
  invisible(ret$Readout)
}
