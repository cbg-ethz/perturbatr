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

#' Select hits from an analyzed data set
#'
#' @export
#' @import data.table
#'
#' @param obj  the data.table to be called
#' @param ...  additional parameters
#'  \itemize{
#'  \item{\emph{hit.ratio} }{ the ratio of siRNAs hits such that a gene counts as hit}
#'  \item{\emph{readout.threshold} }{ the mean readout threshold for all the siRNAs}
#'  \item{\emph{p.value.threshold} }{ p.value threshold for hit selection}
#'  \item{\emph{fdr.threshold} }{ fdr.threshold for hit selection}
#' }
prioritize <- function(obj, ...)
{
  UseMethod("prioritize", obj)
}

#' @export
#' @import data.table
#' @method prioritize svd.analysed.tt
prioritize.svd.analysed.tt <- function(obj, ...)
{
  res <- .select.hits.tt(obj, ...)
  class(res) <- c("svd.prioritized.tt", "svd.prioritized", class(res))
  invisible(res)
}

#' @export
#' @import data.table
#' @method prioritize svd.analysed.hyper
prioritize.svd.analysed.hyper <- function(obj, ...)
{
  res <- .select.hits.hyper(obj, ...)
  class(res) <- c("svd.prioritized.hyper", "svd.prioritized", class(res))
  invisible(res)
}

#' @export
#' @import data.table
#' @method prioritize svd.analysed.pmm
prioritize.svd.analysed.pmm <- function(obj, ...)
{
  pars <- list(...)
  eft  <- ifelse(methods::hasArg(readout.threshold), pars$readout.threshold, .05)
  fdrt <- ifelse(methods::hasArg(fdr.threshold), pars$fdr.threshold, .2)
  res <- .select.hits.pmm(obj, eft, fdrt)
  class(res) <- c("svd.prioritized.pmm", "svd.prioritized", class(res))
  res$fit <- obj
  invisible(res)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr group_by summarize ungroup filter select mutate
#' @importFrom methods hasArg
.select.hits.tt <- function(obj, ...)
{
  # TODO here: what do do with multiple sirnas? same hit criterion as in hyper
  params <- list(...)
  hit.rat <- ifelse(methods::hasArg(hit.ratio), params$hit.ratio, 0.5)
  read.thresh <- ifelse(methods::hasArg(readout.threshold),
                        params$readout.threshold, 0.0)
  p.val.thresh <- ifelse(methods::hasArg(p.value.threshold),
                         params$p.value.threshold, 0.05)
  message(paste("Prioritizing on hit.ratio ", hit.rat,
                ", readout threshold " , read.thresh,
                " and p-value threshold ", p.val.thresh, sep=""))
  res <- dplyr::group_by(obj, Virus, Screen, Library,
                         ScreenType, ReadoutType,
                         Design, Cell,
                         GeneSymbol, Entrez,
                         Plate, RowIdx, ColIdx, siRNAIDs) %>%
    dplyr::mutate(Hit=(Pval <= p.val.thresh & abs(Readout) >= read.thresh)) %>%
    ungroup %>%
    dplyr::group_by(Virus, Screen, Library,
                    ReadoutType, ScreenType,
                    Design, Cell,
                    GeneSymbol, Entrez) %>%
    dplyr::summarize(HitRatio = (sum(Hit)/n()),
                     MeanEffect=mean(Readout),
                     MeanPvalue=mean(Pval)) %>%
    ungroup %>%
    dplyr::filter(HitRatio >= hit.rat)
  res
}

#' @noRd
#' @import data.table
#' @importFrom dplyr group_by summarize ungroup filter select
.select.hits.hyper <- function(obj, ...)
{
  params <- list(...)
  hit.rat      <- ifelse(methods::hasArg(hit.ratio), params$hit.ratio, 0.5)
  message(paste("Prioritizing on hit.ratio ", hit.rat, sep=""))
  res <- dplyr::group_by(obj, Virus, Screen, Library,
                         ScreenType, ReadoutType,
                         Design, Cell,
                         GeneSymbol, Entrez) %>%
    dplyr::summarize(HitRatio        = (base::sum(Hit == TRUE, na.rm=T)/n()),
                     MeanEffect          = base::mean(Readout, na.rm=T),
                     MeanPvalue          = base::mean(Pval)) %>%
    ungroup %>%
    dplyr::filter(HitRatio >= hit.rat)
  res
}

#' @noRd
#' @import data.table
#' @importFrom dplyr filter select group_by mutate summarize full_join
.select.hits.pmm <- function(obj, eft, fdrt)
{
  message(paste0("Prioritizing on fdr.threshold ", fdrt,
                 " and on effect.threshold ", eft ,"."))
  ge <-
    obj$gene.effects %>%
    dplyr::filter(FDR <= fdrt, abs(Effect) >= eft)  %>%
    .[order(-abs(Effect))]
  gpe <- obj$gene.pathogen.effects %>%
    dplyr::filter(FDR <= fdrt, abs(Effect) >= eft)  %>%
    .[order(-abs(Effect))]
  list(gene.hits=ge, gene.pathogen.hits=gpe)
}
