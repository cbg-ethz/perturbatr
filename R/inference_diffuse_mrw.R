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
                do.bootstrap,
                take.largest.component)
{
  message("Diffusion using Markov random walks.")
  diffuse.data <- .do.mrw(hits, adjm, r)
  res          <- diffuse.data$frame

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
             restart.probability    = r,
             delete.nodes.on.degree = delete.nodes.on.degree,
             take.largest.component = take.largest.component),
           dataSet        = hits,
           geneEffects    = res,
           isBootstrapped = is.boot)

  ret
}


#' @noRd
#' @importFrom diffusr random.walk
#' @importFrom assertthat assert_that
.do.mrw <- function(hits, adjm, r)
{
    diffuse.data <- .init.starting.distribution(hits, adjm)
    assertthat::assert_that(
      all(diffuse.data$frame$GeneSymbol == colnames(diffuse.data$adjm)),
      all(diffuse.data$frame$GeneSymbol == rownames(diffuse.data$adjm)))

    mrw          <- diffusr::random.walk(p0=abs(diffuse.data$frame$Effect),
                                         graph=as.matrix(diffuse.data$adjm),
                                         r=r)
    diffuse.data$frame$DiffusionEffect <- mrw
    diffuse.data
}


#' @noRd
#' @import tibble
#' @importFrom dplyr left_join
.init.starting.distribution <- function(hits, adjm)
{
    res <- dplyr::left_join(
      tibble::tibble("GeneSymbol" = colnames(adjm)), hits, by="GeneSymbol")
    res <- dplyr::mutate(
      res, "Effect" = replace(.data$Effect, is.na(.data$Effect), 0))

    list(frame=res, adjm=adjm)
}


#' @noRd
#' @import tibble
#' @importFrom tidyr gather
#' @import foreach
#' @import parallel
#' @import doParallel
#' @importFrom rlang .data
.significance.mrw <- function(bootstrap.hits, adjm, r)
{
    boot.g <- tidyr::gather(bootstrap.hits, "Boot", "Effect",
                            -.data$GeneSymbol)

    doParallel::registerDoParallel(
      ifelse(tolower(Sys.info()['sysname']) %in% c("darwin", "unix"),
             max(1, parallel::detectCores() - 1), 1L))

    lo <- NULL
    li <- foreach::foreach(lo = unique(boot.g$Boot)) %dopar%
    {
      hits  <- dplyr::filter(boot.g, .data$Boot == lo)
      hits  <- dplyr::select(hits, -.data$Boot)
      dd.lo <- .do.mrw(hits, adjm, r)
      dd.lo$frame$boot <- lo
      dd.lo$frame
    }

    doParallel::stopImplicitCluster()

    flat.dat <- dplyr::bind_rows(li)
    flat.dat <- dplyr::select(flat.dat, -.data$Effect)

    ret <- dplyr::group_by(flat.dat, .data$GeneSymbol)
    ret <- ungroup(dplyr::summarise(
      ret, "Mean"= base::mean(.data$DiffusionEffect, na.rm=TRUE)))
    ret <- dplyr::left_join(
      ret,
      tidyr::spread(flat.dat, .data$boot, .data$DiffusionEffect),
      by="GeneSymbol")

    ret
}
