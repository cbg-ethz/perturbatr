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
#' @importFrom tidyr spread
.normalize <- function(obj, normalize, normalize.viability, ...)
{
  res <-
    .normalize.within(obj, normalize, normalize.viability, ...) %>%
    tidyr::spread(ReadoutClass, Readout)
  res
}

#' @noRd
#' @import data.table
#' @importFrom dplyr filter
#' @importFrom methods hasArg
.normalize.within <- function(obj, normalize, normalize.viability, ...)
{
  if (!("Readout" %in% colnames(obj))) stop("No 'Readout' column given!")
  read.dat <- dplyr::filter(obj, ReadoutClass == "Readout")
  via.dat <- dplyr::filter(obj, ReadoutClass == "Viability")
  # pars arguments
  params <- list(...)
  z.score.level <- ifelse(methods::hasArg(z.score.level),
                          params$z.score.level, "plate")
  z.score.ctrl  <- ifelse(methods::hasArg(z.score.ctrl),
                          params$z.score.ctrl, NA_character_)
  poc.ctrl      <- ifelse(methods::hasArg(poc.ctrl), params$poc.ctrl, NA_character_)
  method        <- ifelse(methods::hasArg(method), params$method, "mean")
  background.column <- ifelse(methods::hasArg(background.column),
                              params$background.column,
                              NA)
  background.row <- ifelse(methods::hasArg(methods::background.row),
                           params$background.row, NA)
  # normalize the readout
  message("Normalizing readout")
  read.dat <- .do.normalize.within(obj=read.dat,
                                  normalize=normalize,
                                  method=method,
                                  poc.ctrl=poc.ctrl,
                                  z.score.level=z.score.level,
                                  z.score.ctrl=z.score.ctrl,
                                  background.column=background.column,
                                  background.row=background.row)
  # normalize the viability
  if (normalize.viability)
  {
    message("Normalizing viability")
    via.dat <- .do.normalize.within(obj=via.dat,
                                     normalize=normalize,
                                     method=method,
                                     poc.ctrl=poc.ctrl,
                                     z.score.level=z.score.level,
                                     z.score.ctrl=z.score.ctrl,
                                     background.column=background.column,
                                     background.row=background.row)
  }
  invisible(data.table::rbindlist(list(read.dat, via.dat)))
}

#' @noRd
#' @import data.table
.do.normalize.within <- function(obj, normalize, method, poc.ctrl,
                                 z.score.level, z.score.ctrl,
                                 background.column, background.row)
{
  for (norm in normalize)
  {
    obj <- switch(norm,
                  "log"=.log.norm(obj=obj),
                  "poc"=.poc(obj=obj, method=method, ctrl.gene=poc.ctrl),
                  "z.score"=.z.score(obj=obj, method="default",
                                      level=z.score.level,
                                      ctrl=z.score.ctrl),
                  "robust-z.score"=.z.score(obj=obj, method="robust",
                                            level=z.score.level,
                                            ctrl=z.score.ctrl),
                  "loess"=.loess(obj=obj),
                  "b.score"=.b.score(obj=obj),
                  "background"=.background.correct(
                    obj=obj, background.column=background.column,
                    background.row=background.row, method=method),
                  "qq"=.qq.norm(obj=obj),
                  stop("Give a default normalization method!"))
  }
  obj
}
