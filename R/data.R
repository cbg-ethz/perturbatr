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

#' A sample pan-pathogenic RNAi dataset
#'
#' Example \code{knockout.data} data-set consisting of two-viral screen with
#'  random data.
#'
#' \itemize{
#'     \item Virus  \code{character} names of the viruses
#'     \item Replicate  \code{integer} replicate numer
#'     \item Plate  \code{integer} plate index
#'     \item RowIdx  \code{integer} row index of the well on the plate
#'     \item ColIdx  \code{integer} column index of the well on the plate
#'     \item GeneSymbol  \code{character} HugoSymbol of a gene
#'     \item Entrez  \code{integer} entrez ID
#'     \item ReadoutType  \code{character} readouttype,
#'       such as `GFP` or `Luciferase`
#'     \item Control  \code{integer} coding of controls.
#'       `-1` for negative control, `1` for positive control, `0` for regular
#'        sample
#'     \item Library  \code{character} library type, such as `Ambion`
#'     \item siRNAIDs  \code{character} sirna identifier
#'     \item Screen  \code{character} identifier,
#'       for example `Kinome` or `Genome`
#'     \item Cell  \code{character} cell type, such as `Huh7.5`
#'     \item ScreenType \code{character} screen type, such as `E/R` for
#'        entry/replication
#'     \item Design  \code{character} design of the library, e.g. `pooled`
#'     \item Readout  \code{numeric} readout value
#'     \item ReadoutClass  \code{character} class of the readout, such as
#'       `Readout` or `Viability`
#'     \item NumCells  \code{integer} number of measured cells per well
#'  }
#
#' @author Simon Dirmeier | Computational Biology Group, ETH ZURICH |
#'  \email{simon.dirmeier@@bsse.ethz.ch}
#'
#' @name rnaiscreen
#'
#' @docType data
#' @keywords datasets data
#' @usage data(rnaiscreen)
#' @format A \code{knockout.data} object containing a data.table with 18
#'  columns, each describing a necessary feature
NA
