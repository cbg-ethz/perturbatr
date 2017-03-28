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
#' @importFrom diffusr random.walk
.mrw <- function(hits, bootstrap.hits, adjm, r, graph)
{
  diffuse.data <- .do.mrw(hits, adjm, r)
  if (!is.null(bootstrap.hits))
  {
    pvals <- .significance.mrw(bootstrap.hits, adjm, r)
    res   <- dplyr::left_join(diffuse.data$frame, pvals, by="GeneSymbol")
  }
  else
  {
    res  <- diffuse.data$frame
  }
  li  <- list(diffusion=dplyr::select(res, GeneSymbol, Effect, DiffusionEffect),
              diffusion.model=res,
              lmm.hits=hits,
              graph=graph)
  li
}

#' @noRd
#' @importFrom diffusr random.walk
.do.mrw <- function(hits, adjm, r)
{
  diffuse.data <- .init.starting.distribution(hits, adjm)
  mrw <- diffusr::random.walk(abs(diffuse.data$frame$Effect),
                              as.matrix(diffuse.data$adjm), r)
  diffuse.data$frame$DiffusionEffect <- mrw
  diffuse.data
}

#' @noRd
#' @import data.table
#' @importFrom dplyr left_join
.init.starting.distribution <- function(hits, adjm)
{
  res  <- dplyr::left_join(data.table::data.table(GeneSymbol=colnames(adjm)),
                           hits, by="GeneSymbol")
  data.table::setDT(res)[is.na(Effect), Effect := 0]
  list(frame=res, adjm=adjm)
}

#' @noRd
#' @importFrom tidyr gather
#' @import foreach
#' @import parallel
#' @import doParallel
.significance.mrw <- function(bootstrap.hits, adjm, r)
{
  boot.g <- tidyr::gather(bootstrap.hits, Boot, Effect, -GeneSymbol) %>%
    as.data.table
  li <- list()
  # TODO: parallel
  foreach::foreach(lo=unique(boot.g$Boot)) %do%
  {
    hits <- dplyr::filter(boot.g, Boot==lo)
    hits <- dplyr::select(hits, -Boot)
    dd.lo <- .do.mrw(hits, adjm, r)
    dd.lo$frame$boot <- lo
    li[[lo]] <- dd.lo$frame
  }
  # TODO: how to do hypothesis test here?
  # -> calculate means -> calculate alphas -> calculate variance -> calculate alpha_0
  # prolly as with lmm, however ttest is a ad choice here, since it depends on
  # the graph size
  flat.dat <- do.call("rbind", lapply(li, function(e) e)) %>%
    dplyr::select(-Effect)
  dat <- flat.dat  %>%
    dplyr::group_by(GeneSymbol) %>%
    dplyr::summarise(MeanBootstrapDiffusionEffect=
                       mean(DiffusionEffect, na.rm=T)) %>%
    ungroup
  dat <- dplyr::left_join(dat,
                          tidyr::spread(flat.dat, boot, DiffusionEffect),
                          by="GeneSymbol")
  dat
}
