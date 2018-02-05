# perturbR: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2018 Simon Dirmeier
#
# This file is part of perturbR
#
# perturbR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# perturbR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with perturbR. If not, see <http://www.gnu.org/licenses/>.


#' @noRd
#' @import data.table
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
.z.score <- function(obj,
                     method = c("default", "robust"),
                     level  = c("plate", "control"),
                     ctrl)
{
  method <- match.arg(method)
  level  <- match.arg(level)
  message(paste("Calculating", method, "z-scores on", level))
  if (!is.na(ctrl)) message(paste("...normalizing on", ctrl))
  z.score.data <-
    dplyr::group_by(obj, Condition, Screen, Library,
                    ReadoutType, ScreenType, ReadoutClass,
                    Design, Cell, Replicate, Plate)
  if (level == "plate")
  {
    z.score.data <-
      dplyr::mutate(z.score.data, Readout = .z.score.plate(Readout, method))
  }
  else if (level == "control")
  {
    check <- dplyr::filter(z.score.data, Control == -1) %>%
      dplyr::mutate(n=n())
    if (length(unique(check$n))  > 1)
      warning(paste("Found unequal number of negative controls on plates:",
                    paste(unique(check$n), collapse=", ")))
    z.score.data <-
      dplyr::mutate(z.score.data,
                    Readout = .z.score.control(Readout, method,
                                               Control, GeneSymbol, ctrl))
  }
  dplyr::ungroup(z.score.data)
  invisible(z.score.data)
}

#' @noRd
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom stats mad
.z.score.plate <- function(obj, method)
{
  val <- switch(method,
                "default"=((obj - mean(obj, na.rm=TRUE)) /
                             stats::sd(obj, na.rm=TRUE)),
                "robust" =((obj - stats::median(obj, na.rm=TRUE)) /
                             stats::mad(obj, na.rm=TRUE)),
                stop("No correct method given!"))
  invisible(val)
}

#' @noRd
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom stats mad
.z.score.control <- function(re, method, control, genes, ctrl.gene)
{
  genes <- tolower(genes)
  ctrl.gene <- tolower(ctrl.gene)
  cont.idx <- which(control == -1)
  if (!is.na(ctrl.gene)) {
    cont.idx <- which(control == -1 & genes == ctrl.gene)
    if (length(cont.idx) == 0) stop("No controls found for criteria!")
  }
  if (is.na(cont.idx[1]))
  {
    warning("There are no negative controls on your plate!
            Stopping control-wise normalization and using plate-wise instead.")
    val <- .z.score.plate(re, method)
  }
  else
  {
    if (length(cont.idx) < 3)
      warning(paste("You are normalizing with z-scores and only using",
                    length(cont.idx),"values!"))
    val <- switch(method,
                  "default"=((re - mean(re[cont.idx], na.rm=TRUE)) /
                               stats::sd(re[cont.idx], na.rm=TRUE)),
                  "robust" =((re - stats::median(re[cont.idx], na.rm=TRUE)) /
                               stats::mad(re[cont.idx], na.rm=TRUE)),
                  stop("No correct method given!"))
  }
  val
}
