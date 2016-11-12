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
.z.score <- function(obj, method=c("default", "robust"),
                     level=c("plate", "control"), ctrl)
{
  method <- match.arg(method)
  level  <- match.arg(level)
  message(paste("Calculating z-scores on", level))
  if (!is.na(ctrl)) message(paste("...normalizing on", ctrl))
  z.score.data <-
    dplyr::group_by(obj, Virus, Screen, Library,
                    ReadoutType, InfectionType, ReadoutClass,
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
                "default"=((obj - mean(obj, na.rm=T)) /
                             stats::sd(obj, na.rm=T)),
                "robust" =((obj - stats::median(obj, na.rm=T)) /
                             stats::mad(obj, na.rm=T)),
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
                  "default"=((re - mean(re[cont.idx], na.rm=T)) /
                               stats::sd(re[cont.idx], na.rm=T)),
                  "robust" =((re - stats::median(re[cont.idx], na.rm=T)) /
                               stats::mad(re[cont.idx], na.rm=T)),
                  stop("No correct method given!"))
  }
  val
}
