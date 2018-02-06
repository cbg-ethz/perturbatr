# perturbR: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2018 Simon Dirmeier
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
# along with perturbR. If not, see <http://www.gnu.org/licenses/>.


#' @include class_data.R
#' @include class_analysed.R


#' @title Plot a \code{perturbation.mrw.diffusion.analysed} object
#'
#' @description Plots a \code{perturbation.mrw.diffusion.analysed} object by
#' iteratively expanding the neighbors of the hits form the diffusion
#'
#' @method plot perturbation.diffusion.analysed
#' @export
#'
#' @import data.table
#' @import grid
#' @importFrom gridExtra tableGrob ttheme_minimal grid.table
#'
#' @param x  a \code{perturbation.diffusion.analysed} object
#' @param cnt  number of genes to be shown
#' @param ...  additional params
#'
#' @return returns a plot object
#'
plot.perturbation.diffusion.analysed <- function(x, cnt=10, ...)
{
  tt <- gridExtra::ttheme_minimal(
    base_family="Helvetica",
    core=list(bg_params = list(fill="white")))

  x@.data %>%
    .[order(-DiffusionEffect)] %>%
    .[1:cnt] %>%
    as.data.frame %>%
    gridExtra::grid.table(theme=tt)
}
