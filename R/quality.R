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

#' Calculates per-plate, per-replicate and screen quality scores
#'
#' @export
#' @import data.table
#'
#' @param obj  the object for which quality scores are calculates
#' @param ...  additional parameters
quality <- function(obj, ...) UseMethod("quality")


#' @export
#' @noRd
#' @import data.table
#' @importFrom dplyr filter
quality.svd.raw <-
function
(
  obj,
  ...
)
{
  ret <- obj %>%
    dplyr::filter(ReadoutClass == "Readout")
  quality.svd.data(ret)
}

#' @export
#' @noRd
#' @import data.table
#' @importFrom dplyr select
quality.svd.data <-
function
(
  obj,
  ...
)
{
  q <- .quality(obj)
  res <- dplyr::select(obj, Virus, Screen, Library, Cell,
                       InfectionType, ReadoutType,
                       Replicate, Plate, Control, Readout)
  res <- data.table::as.data.table(res)
  ret <- list(quality=q, data=res)
  class(ret) <- "svd.quality"
  invisible(ret)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr group_by summarize
.quality <-
function
(
    obj
)
{
  q.plate <- obj %>%
    dplyr::group_by(Virus, Screen, Library, Replicate,
                    Plate, ReadoutType, InfectionType, Cell, Design) %>%
    dplyr::summarize(z.fac.control=.z.factor(Readout, Control),
                     ssmd=.ssmd(Readout, Control),
                     z.fac.plate=.z.factor(Readout, Control, "plate")) %>%
    ungroup
  q.rep <- obj %>%
    dplyr::group_by(Virus, Screen, Library, Replicate,
                    ReadoutType, InfectionType, Cell, Design) %>%
    dplyr::summarize(z.fac.control=.z.factor(Readout, Control),
                     ssmd=.ssmd(Readout, Control),
                     z.fac.plate=.z.factor(Readout, Control, "plate")) %>%
    ungroup
  q.screen <- obj %>%
    dplyr::group_by(Virus, Screen, Library,
                    ReadoutType, InfectionType, Cell, Design) %>%
    dplyr::summarize(z.fac.control=.z.factor(Readout, Control),
                     ssmd=.ssmd(Readout, Control),
                     z.fac.plate=.z.factor(Readout, Control, "plate")) %>%
    ungroup

  list(plate.quality=q.plate,
       replicate.quality=q.rep,
       screen.quality=q.screen)
}

#' @noRd
#' @importFrom stats sd
.z.factor <- function(read, ctrl, level=c("control", "plate"))
{
    level <- match.arg(level)
    neg.crtl <- which(ctrl == -1)
    if (level == "control")
    {
      ot.crtl <- which(ctrl == 1)
    }
    else if (level=="plate")
    {
      ot.crtl <- which(ctrl == 0)
    }
    else
    {
      stop("Please provide a standard method")
    }
    z.fac <- -Inf
    if (length(ot.crtl) == 0| length(neg.crtl) == 0)
    {
      warning("Could not find enough control indexes. Returning neg. infinity!")
    }
    else
    {
      den <- 1 - ( 3 * stats::sd(read[ot.crtl]) +
                     3 *  stats::sd(read[neg.crtl]))
      nom <- abs(base::mean(read[ot.crtl]) - base::mean(read[neg.crtl]))
      z.fac <-  den/nom
    }
    z.fac
}

#' @noRd
#' @importFrom stats var
.ssmd <- function(read, ctrl)
{
  pos.crtl <- which(ctrl == 1)
  neg.crtl <- which(ctrl == -1)
  ssmd <- -Inf
  if (length(pos.crtl) == 0| length(neg.crtl) == 0)
  {
    warning("Could not find enough control indexes. Returning neg. infinity!")
  }
  else
  {
    pos <- read[pos.crtl]
    neg <- read[neg.crtl]
    ssmd <- abs(base::mean(pos) - base::mean(neg)) /
      (stats::var(pos) + stats::var(neg))
  }
  ssmd
}

