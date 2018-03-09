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
# along with perturbatr. If not, see <http://www.gnu.org/licenses/>.


#' @include class_data.R
#' @include class_analysed.R


#' @title Plot a \code{NetworkAnalysedPerturbationData} object
#'
#' @description Creates a table of the gene ranking of a
#' \code{NetworkAnalysedPerturbationData} object
#'
#' @method plot NetworkAnalysedPerturbationData
#' @export
#'
#' @import data.table
#' @import grid
#' @importFrom gridExtra tableGrob ttheme_minimal grid.table
#'
#' @param x  a \code{NetworkAnalysedPerturbationData} object
#' @param cnt  number of genes to be shown
#' @param ...  additional parameters
#'
#' @return returns a table if the first \code{cnt} highest ranked genes
#'
plot.NetworkAnalysedPerturbationData <- function(x, cnt=10, ...)
{
  tt <- gridExtra::ttheme_minimal(
    base_family="Helvetica",
    core=list(bg_params = list(fill="white")))

  geneEffects(x) %>%
    .[order(-DiffusionEffect)] %>%
    .[seq(cnt)] %>%
    as.data.frame() %>%
    gridExtra::grid.table(theme=tt)
}
