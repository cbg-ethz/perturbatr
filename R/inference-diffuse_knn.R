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
#' @import data.table igraph
#' @importFrom dplyr filter mutate
.knn <- function(hits, bootstrap.hits, adjm, node.start.count,
                 search.depth, graph)
{
  # TODO diff on bootstrap hits
  diffuse.data <- .init.starting.indexes(hits, adjm, node.start.count)

  neighs <- diffusr::nearest.neighbors(
    as.integer(dplyr::filter(diffuse.data$frame, Select==T)$Idx),
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

  genes.found <- dplyr::select(flat.dat, Start, GeneSymbol) %>%
    unlist %>% unname %>% unique

  genes.f.table <- data.table::data.table(GeneSymbol=genes.found) %>%
    dplyr::left_join(res, by="GeneSymbol")

  li <- list(data=res,
             neighbors=flat.dat,
             graph=graph,
             genes.found=genes.f.table)

  li
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
