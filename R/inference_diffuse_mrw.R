# perturbatr: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2018 Simon Dirmeier
#
# This file is part of perturbatr
#
# perturbatr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# perturbatr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with perturbatr If not, see <http://www.gnu.org/licenses/>.


#' @include class_data.R
#' @include class_analysed.R


#' @noRd
#' @import tibble
#' @importFrom diffusr random.walk
#' @importFrom methods new
mrw <- function(hits,
                mod,
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


  is.boot <- if (!is.null(bootstrap.hits) && do.bootstrap)
  {
      boot.intrvls <- .significance.mrw(bootstrap.hits, adjm, r)
      res          <- dplyr::left_join(diffuse.data$frame,
                                       boot.intrvls,
                                       by="GeneSymbol")
      TRUE
  } else FALSE

  ret <- methods::new("NetworkAnalysedPerturbationData",
           graph           = graph,
           params          = list(
             restart.probability     = r,
             delete.nodes.on.degree = delete.nodes.on.degree),
           dataSet        = hits,
           geneEffects    = res,
           isBootstrapped = is.boot)

  ret
}


#' @noRd
#' @importFrom diffusr random.walk
.do.mrw <- function(hits, adjm, r)
{
    diffuse.data <- .init.starting.distribution(hits, adjm)
    mrw          <- diffusr::random.walk(abs(diffuse.data$frame$Effect),
                                         as.matrix(diffuse.data$adjm),
                                         r)
    diffuse.data$frame$DiffusionEffect <- mrw
    diffuse.data
}


#' @noRd
#' @import tibble
#' @importFrom dplyr left_join
.init.starting.distribution <- function(hits, adjm)
{
    res <- dplyr::left_join(
      tibble::tibble(GeneSymbol=colnames(adjm)), hits, by="GeneSymbol")
    res <- dplyr::mutate(res, Effect = replace(is.na(Effect), 0, Effect))

    list(frame=res, adjm=adjm)
}


#' @noRd
#' @import tibble
#' @importFrom tidyr gather
#' @import foreach
#' @import parallel
#' @import doParallel
.significance.mrw <- function(bootstrap.hits, adjm, r)
{
    boot.g <- tidyr::gather(bootstrap.hits, Boot, Effect, -GeneSymbol)

    doParallel::registerDoParallel(
      ifelse(tolower(Sys.info()['sysname']) %in% c("darwin", "unix"),
             max(1, parallel::detectCores() - 1), 1L))

    li <- foreach::foreach(lo=unique(boot.g$Boot)) %dopar%
    {
      hits  <- dplyr::filter(boot.g, Boot == lo)
      hits  <- dplyr::select(hits, -Boot)
      dd.lo <- .do.mrw(hits, adjm, r)
      dd.lo$frame$boot <- lo
      dd.lo$frame
    }

    doParallel::stopImplicitCluster()

    flat.dat <- dplyr::bind_rows(li) %>%
      tibble::as.tibble() %>%
      dplyr::select(-Effect)
    ret <- flat.dat  %>%
      dplyr::group_by(GeneSymbol) %>%
      dplyr::summarise(Mean=mean(DiffusionEffect, na.rm=TRUE)) %>%
      ungroup() %>%
      dplyr::left_join(tidyr::spread(flat.dat, boot, DiffusionEffect),
                       by="GeneSymbol")

    ret
}
