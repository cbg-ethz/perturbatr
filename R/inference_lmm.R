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


#' @title Fit an LMM to the data and calculate local false discovery rates.
#'
#' @description \code{lmm} TODO
#'
#' @export
#' @docType methods
#' @rdname lmm-methods
#'
#' @import data.table
#'
#' @param obj  an svd.data object
#' @param drop  boolean flag if all entries should be dropped that are not
#'  found in every virus
#' @param weights a list of weights
#' @param rel.mat.path  the (optional) path to a target relation matrix that is
#'  going to be used for
#' @param bootstrap.cnt  the number of loocv runs you want to do in order to
#'  estimate a significance level for the gene effects
#' @param ignore  ignore siRNAs that have been seen only once per group
#' @param ...  additional parameters
setGeneric(
  "lmm",
  function(obj,
           drop=T,
           weights=NULL,
           rel.mat.path=NULL,
           bootstrap.cnt=F,
           ignore=1, ...)
  {
    standardGeneric("lmm")
  },
  package="knockout"
)

#' @rdname lmm-methods
#' @aliases lmm,knockout.lmm.data-method
#' @import data.table
setMethod(
  "lmm",
  signature = signature(obj="knockout.data"),
  function(obj,
           drop=T,
           weights=NULL,
           rel.mat.path=NULL,
           bootstrap.cnt=F,
           ignore=1,
           ...)
  {
    .check.data(obj)
    res <- .lmm.knockout.data(obj@.data, bootstrap.cnt)
    ret <- new("knockout.analysed",
               .inference=.inference.types()$MIXED.MODEL,
               .data=res)
    ret
  }
)

#' @rdname lmm-methods
#' @aliases lmm,knockout.lmm.data-method
#' @import data.table
setMethod(
  "lmm",
  signature = signature(obj="knockout.lmm.data"),
  function(obj,
           drop=T,
           weights=NULL,
           rel.mat.path=NULL,
           bootstrap.cnt=F,
           ignore=1,
           ...)
  {
    res <- .lmm.model.data(obj@.data, bootstrap.cnt)
    ret <- new("knockout.analysed",
               .inference=.inference.types()$MIXED.MODEL,
               .data=res)
    ret
  }
)

#' @noRd
#' @import data.table
.lmm.knockout.data <- function(obj, drop,
                          weights=NULL, rel.mat.path=NULL,
                          bootstrap.cnt, ignore, ...)
{
  # init the data table for the LMM
  md <- .set.lmm.matrix(obj, drop, ignore, weights, rel.mat.path)
  .lmm.model.data(md, bootstrap.cnt)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr mutate full_join select group_by
.lmm.model.data <- function(md, bootstrap.cnt)
{
  if (is.numeric(bootstrap.cnt) & bootstrap.cnt < 10 & bootstrap.cnt >= 1)
  {
    stop("Please use at least 10 bootstrap runs (better 100/1000).")
  }
  # save gene control mappings
  gene.control.map <-
    dplyr::select(md, GeneSymbol, Control) %>%
    unique %>%
    dplyr::mutate(GeneSymbol=as.character(GeneSymbol))
  mult.gen.cnt <- (gene.control.map %>%
                     dplyr::group_by(GeneSymbol) %>%
                     dplyr::mutate(cnt=n()))$cnt %>%
    unique
  if (length(mult.gen.cnt) != 1)
  {
    warning("Found multiple gene-control entries for several genes,
            i.e. several genes are both control and not control!")
  }
  # fit the LMM
  message("Fitting LMM")
  fit.lmm <- .lmm(md)
  ref <- .ranef(fit.lmm)
  # calculate fdrs
  gp.fdrs <- .fdr(ref$gene.pathogen.effects)
  # TODO make this consistent to previous FDR
  if (is.numeric(bootstrap.cnt) & bootstrap.cnt >= 10)
  {
    message("Bootstrap for significance estimation")
    ge.fdrs <- .lmm.significant.hits(md, bootstrap.cnt)
  }
  else ge.fdrs <- data.table(GeneSymbol=ref$gene.effects$GeneSymbol, FDR=NA_real_)
  # set together the gene/fdr/effects and the mappings
  gene.effects <- dplyr::full_join(ref$gene.effects, gene.control.map,
                                   by="GeneSymbol") %>%
    dplyr::full_join(dplyr::select(ge.fdrs, GeneSymbol, FDR), by="GeneSymbol")
  gene.path.effs <- dplyr::full_join(gp.fdrs$gene.pathogen.matrix,
                                     gene.control.map, by="GeneSymbol")
  ret <- list(gene.effects=gene.effects,
              gene.pathogen.effects=gene.path.effs,
              infection.effects=ref$infection.effects,
              model.data=md, fit=list(model=fit.lmm,
                                      gene.pathogen.fdrs=gp.fdrs$fdrs,
                                      gene.fdrs=ge.fdrs))
  ret
}

#' @noRd
.init.formula <- function()
{
  frm.str <- paste0("Readout ~ Virus + ",
                    "(1 | GeneSymbol) + (1 | Virus:GeneSymbol) + ",
                    "(1 | ScreenType) + (1 | Virus:ScreenType)")
  frm.str
}

#' @noRd
#' @import data.table
#' @importFrom lme4 lmer
#' @importFrom stats as.formula
.lmm <- function(md)
{
  lme4::lmer(stats::as.formula(.init.formula()),
             data = md, weights = md$Weight, verbose = F)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr mutate select
#' @importFrom lme4 ranef
.ranef <-  function(fit.lmm)
{
  random.effects <- lme4::ranef(fit.lmm)
  # create the data table with gene effects
  ge <- data.table::data.table(
    Effect = random.effects[["GeneSymbol"]][,1],
    GeneSymbol = as.character(rownames(random.effects[["GeneSymbol"]])))
  # create the data.table with gene-pathogen effects
  gpe <- data.table::data.table(
    gpe = random.effects[["Virus:GeneSymbol"]][,1],
    GenePathID =
      as.character(rownames(random.effects[["Virus:GeneSymbol"]]))) %>%
    dplyr::mutate(GeneSymbol = sub("^.+:", "", GenePathID))
  # table for infection types
  ie <- data.table::data.table(
    Effect = random.effects[["ScreenType"]][,1],
    ScreenType = as.character(rownames(random.effects[["ScreenType"]])))
  # table for virus-infection types
  # ipe <- data.table::data.table(
  #   ie = random.effects[["Virus:ScreenType"]][,1],
  #   InfectionPathID = as.character(rownames(random.effects[["Virus:ScreenType"]])))
  # create the table with gene-pathogen effects
  ga <- base::merge(gpe, ge, by = "GeneSymbol") %>%
    dplyr::mutate(Virus = sub(":.+$", "", GenePathID),
                  GeneVirusEffect = Effect + gpe) %>%
    dplyr::select(-GenePathID, -gpe, -Effect)
  list(gene.effects=ge,
       gene.pathogen.effects=ga,
       infection.effects=ie)
}
