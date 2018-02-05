# perturbR: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2015 - 2016 Simon Dirmeier
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
# along with perturbR If not, see <http://www.gnu.org/licenses/>.



#' Calculates per-plate, per-replicate and screen quality scores
#'
#' @export
#' @import data.table
#'
#' @param obj  the object for which quality scores are calculates
#' @param ...  additional parameters
#'
#' @return returns a \code{perturbation.quality} object
quality <- function(obj, ...)
{
  UseMethod("quality")
}

#' @export
#' @method quality perturbation.raw.data
#' @import data.table
#' @importFrom dplyr filter
quality.perturbation.raw.data <-function(obj, ...)
{
  obj@.data <- obj@.data %>%
    dplyr::filter(ReadoutClass == "Readout")
  quality.perturbation.normalized.data(obj)
}

#' @export
#' @method quality perturbation.normalized.data
#' @import data.table
#' @importFrom dplyr select
quality.perturbation.normalized.data <- function(obj, ...)
{
  q <- .quality(obj@.data)
  new("perturbation.quality", .quality=q, .data=obj@.data)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr group_by summarize
.quality <- function(obj)
{
  q.plate <- obj %>%
    dplyr::group_by(Virus, Screen, Library, Replicate,
                    Plate, ReadoutType, ScreenType, Cell, Design) %>%
    dplyr::summarize(z.fac.control = .z.factor(Readout, Control),
                     ssmd          = .ssmd(Readout, Control),
                     z.fac.plate   = .z.factor(Readout, Control, "plate")) %>%
    ungroup
  q.rep <- obj %>%
    dplyr::group_by(Virus, Screen, Library, Replicate,
                    ReadoutType, ScreenType, Cell, Design) %>%
    dplyr::summarize(z.fac.control=.z.factor(Readout, Control),
                     ssmd=.ssmd(Readout, Control),
                     z.fac.plate=.z.factor(Readout, Control, "plate")) %>%
    ungroup
  q.screen <- obj %>%
    dplyr::group_by(Virus, Screen, Library,
                    ReadoutType, ScreenType, Cell, Design) %>%
    dplyr::summarize(z.fac.control=.z.factor(Readout, Control),
                     ssmd=.ssmd(Readout, Control),
                     z.fac.plate=.z.factor(Readout, Control, "plate")) %>%
    ungroup

  list(plate.quality     = q.plate,
       replicate.quality = q.rep,
       screen.quality    = q.screen)
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
      den <- 1 -
        ( 3 * stats::sd(read[ot.crtl]) + 3 *  stats::sd(read[neg.crtl]))
      nom <- abs(mean(read[ot.crtl]) - mean(read[neg.crtl]))
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
    ssmd <- abs(mean(pos) - mean(neg)) / (stats::var(pos) + stats::var(neg))
  }
  ssmd
}
