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
.b.score <- function(obj)
{
  message("Calculating b-scores!")
  b.score.data <-
    dplyr::group_by(obj, Virus, Screen, Library,
                    ReadoutType, ScreenType, ReadoutClass,
                    Design, Cell,
                    Replicate, Plate) %>%
    dplyr::mutate(Readout = .bscore.plate(RowIdx, ColIdx, Readout)) %>%
    ungroup
  invisible(b.score.data)
}

#' @noRd
#' @importFrom stats medpolish
#' @importFrom stats median
#' @importFrom stats mad
.bscore.plate <- function(rows, cols, readouts)
{
  nrow <- max(rows)
  ncol <- max(cols)
  fr <- data.frame(rows, cols, readouts, ord=1:length(readouts))
  fr <- fr[order(fr$rows, fr$cols),]
  # this is needed in case the data is not sorted
  readout.mat <- matrix(fr$readouts, nrow, ncol, byrow=TRUE)
  res <- rep(NA_real_, length(readouts))
  tryCatch ({
    med.pol <- stats::medpolish(readout.mat, maxiter=100,
                                na.rm=TRUE, trace.iter=FALSE)$residuals
    res <- as.vector(t(med.pol)) / stats::mad(med.pol, na.rm=TRUE)
  }, warning = function(war) { message(paste("WARNING: ", war)); },
  error = function(err)   { message(paste("ERROR: ", err));   }
  )
  if (all(is.na(res))) res <- fr$readouts
  fr$res <- res
  fr <- fr[order(fr$ord),]
  invisible(fr$res)
}
