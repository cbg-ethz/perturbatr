# knockdown: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2015 - 2016 Simon Dirmeier
#
# This file is part of knockdown
#
# knockdown is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# knockdown is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with knockdown. If not, see <http://www.gnu.org/licenses/>.


#' @include class_knockdown_data.R
#' @include class_knockdown_analysed.R


#' @noRd
#' @import data.table igraph
#' @importFrom dplyr filter mutate
knn <- function(hits,
                mod,
                delete.nodes.on.degree,
                bootstrap.hits,
                adjm,
                node.start.count,
                search.depth,
                graph,
                do.bootstrap)
{

  message("Diffusion using kNN.")

  # TODO diff on bootstrap hits
  # TODO diff on bootstrap hits
  # TODO diff on bootstrap hits

  diffuse.data <- .init.starting.indexes(hits, adjm, node.start.count)
  indexes <- dplyr::filter(diffuse.data$frame, Select == TRUE) %>%
    dplyr::select(Idx) %>%
    unlist %>%
    unname %>%
    as.integer

  neighs <- diffusr::nearest.neighbors(
    as.integer(indexes),
    as.matrix(diffuse.data$adjm), search.depth)

  res <- diffuse.data$frame

  names(neighs) <-
    (diffuse.data$frame %>%
       dplyr::filter(Idx %in% as.integer(names(neighs))))$GeneSymbol

  # Add the effects of the neighbors
  neighs <- lapply(neighs, function(e) {
    dplyr::left_join(data.table(Idx=e), res, by="Idx") %>%
      dplyr::select(GeneSymbol, Effect)}
  )

  flat.dat <- do.call(
    "rbind",
    lapply(1:length(neighs),
           function(e) data.table(Start=names(neighs)[e], neighs[[e]]))
  )

  # genes.found <- dplyr::select(flat.dat, Start, GeneSymbol) %>%
  #   unlist %>% unname %>% unique
  #
  # genes.f.table <- data.table::data.table(GeneSymbol=genes.found) %>%
  #   dplyr::left_join(res, by="GeneSymbol")


  ret <- new("knockdown.knn.diffusion.analysed",
             .graph           = graph,
             .neighbors       = data.table::as.data.table(neighs),
             .initial.model   = mod,
             .params          = list(
               node.start.count = node.start.count,
               search.depth    = search.depth,
               delete.nodes.on.degree = delete.nodes.on.degree),
             .inference       = .inference.types()$NN.DIFFUSION,
             .data            = data.table::as.data.table(res),
             .is.bootstrapped = FALSE)

  ret
}

#' @noRd
#' @import data.table
#' @importFrom dplyr left_join mutate select
.init.starting.indexes <- function(hits, adjm, node.start.count)
{
  res  <- dplyr::left_join(data.table::data.table(GeneSymbol=colnames(adjm)),
                           hits, by="GeneSymbol")
  data.table::setDT(res)[is.na(Effect), Effect := 0]
  res <- res %>%
    .[order(-abs(Effect))] %>%
    dplyr::mutate(Idx = 1:.N) %>%
    dplyr::mutate(Select = (Idx <= node.start.count))

  data.table::setDT(res)[Effect == 0  , Select := FALSE]
  data.table::setDT(res)[is.na(Select), Select := FALSE]

  list(frame=res, adjm=adjm)
}
