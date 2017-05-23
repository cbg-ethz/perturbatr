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


#' @include class_knockout_data.R
#' @include class_knockout_analysed.R


#' @noRd
#' @import data.table
#' @importFrom diffusr random.walk
mrw <- function(hits,
                mod,
                do.bootstrap,
                bootstrap.hits,
                delete.nodes.on.degree,
                r,
                adjm,
                graph,
                do.bootstrap)
{
  message("Diffusion using Markov random walks.")
  diffuse.data <- .do.mrw(hits, adjm, r)
  res  <- diffuse.data$frame

  if (!is.null(bootstrap.hits) && do.bootstrap)
  {
    boot.intrvls <- .significance.mrw(bootstrap.hits, adjm, r)
    res   <- dplyr::left_join(diffuse.data$frame,
                              boot.intrvls,
                              by="GeneSymbol")
  }

  ret <- new("knockout.diffusion.analysed",
             .graph           = graph,
             .initial.model   = mod,
             .parameters      = list(
               restart.probaility     = r,
               delete.nodes.on.degree = delete.nodes.on.degree),
             .inference       = .inference.types()$MRW.DIFFUSION,
             .data            = data.table::as.data.table(res),
             .is.bootstrapped = mod@.is.bootstrapped
  )

  ret
}

#' @noRd
#' @importFrom diffusr random.walk
.do.mrw <- function(hits, adjm, r)
{
  diffuse.data <- .init.starting.distribution(hits, adjm)
  mrw          <- diffusr::random.walk(
    abs(diffuse.data$frame$Effect),
    as.matrix(diffuse.data$adjm),
    r)
  diffuse.data$frame$DiffusionEffect <- mrw
  diffuse.data
}

#' @noRd
#' @import data.table
#' @importFrom dplyr left_join
.init.starting.distribution <- function(hits, adjm)
{
  res  <- dplyr::left_join(data.table::data.table(GeneSymbol=colnames(adjm)),
                           hits,
                           by="GeneSymbol")
  data.table::setDT(res)[is.na(Effect), Effect := 0]
  list(frame=res, adjm=adjm)
}

#' @noRd
#' @import data.table
#' @importFrom tidyr gather
#' @import foreach
#' @import parallel
#' @import doParallel
.significance.mrw <- function(bootstrap.hits, adjm, r)
{
  boot.g <- data.table::as.data.table(
    tidyr::gather(bootstrap.hits, Boot, Effect, -GeneSymbol))

  # init multicore
  doParallel::registerDoParallel(
    ifelse(tolower(Sys.info()['sysname']) %in% c("darwin", "unix"),
           max(1, parallel::detectCores() - 1), 1L)
  )
  # to diffusion
  li <- foreach::foreach(lo=unique(boot.g$Boot)) %dopar%
  {
    hits  <- dplyr::filter(boot.g, Boot==lo)
    hits  <- dplyr::select(hits, -Boot)
    dd.lo <- .do.mrw(hits, adjm, r)
    dd.lo$frame$boot <- lo
    dd.lo$frame
  }
  # stop multiprocessing
  doParallel::stopImplicitCluster()

  # TODO: how to do hypothesis test here?
  # prolly using some dirichlet-kind of thing
  # -> calculate means -> calculate alphas ->
  # -> calculate variance -> calculate alpha_0
  flat.dat <- data.table::rbindlist(li) %>% dplyr::select(-Effect)
  ret <- flat.dat  %>%
    dplyr::group_by(GeneSymbol) %>%
    dplyr::summarise(Mean=mean(DiffusionEffect, na.rm=T)) %>%
    ungroup %>%
    dplyr::left_join(tidyr::spread(flat.dat, boot, DiffusionEffect),
                     by="GeneSymbol")

  ret
}
