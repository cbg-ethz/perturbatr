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
# along with perturbR If not, see <http://www.gnu.org/licenses/>.


#' @include class_data.R


#' @title Jointly analyse multiple genetic perturbation screens
#'
#' @description Analyse multiple different genetic perturbation screens at once
#'  using a hierarchical model. The model estimates general relative effect
#'  sizes for genes across all experiments. This could for instance be a
#'  pan-pathogenic host factor, i.e. a gene that decisively impacts the
#'  life-cycle of multiple pathogens.
#'
#' @export
#' @docType methods
#' @rdname hm-methods
#'
#' @import data.table
#'
#' @param obj  an svd.data object
#' @param drop  boolean flag if all entries should be dropped that are not
#'  found in every Condition
#' @param weights a list of weights
#' @param rel.mat.path  the (optional) path to a target relation matrix that is
#'  going to be used for
#' @param bootstrap.cnt  the number of bootstrap runs you want to do in order
#'  to estimate a significance level for the gene effects
#' @param ignore  ignore siRNAs that have been seen only once per group
#' @param effect.size  the relative effect size used for hit prioritization
#' @param qval.threshold  the q-value threshold used for hit prioritization
#'  if bootstrap.cnt is set
#'
#' @return returns a \code{perturbation.hm.analysed} object
#'
setGeneric(
  "hm",
  function(obj,
           drop=TRUE,
           weights=NULL,
           rel.mat.path=NULL,
           bootstrap.cnt=0,
           ignore=1,
           effect.size=0.05,
           qval.threshold=.2)
  {
    standardGeneric("hm")
  },
  package="perturbation"
)

#' @rdname hm-methods
#' @aliases hm,perturbation.normalized.data-method
#' @import data.table
setMethod(
  "hm",
  signature = signature(obj="perturbation.normalized.data"),
  function(obj,
           drop=TRUE,
           weights=NULL,
           rel.mat.path=NULL,
           bootstrap.cnt=0,
           ignore=1,
           effect.size=0.05,
           qval.threshold=.2)
  {
    md <- set.hm.model.data(obj, drop, ignore, weights, rel.mat.path)
    hm(md, drop, weights, rel.mat.path, bootstrap.cnt, ignore,
        effect.size, qval.threshold)
  }
)

#' @rdname hm-methods
#' @aliases hm,perturbation.hm.data-method
#' @import data.table
setMethod(
  "hm",
  signature = signature(obj="perturbation.hm.data"),
  function(obj,
           drop=TRUE,
           weights=NULL,
           rel.mat.path=NULL,
           bootstrap.cnt=0,
           ignore=1,
           effect.size=0.01,
           qval.threshold=.2)
  {
    res     <- .hm.model.data(obj, bootstrap.cnt)
    priorit <- .prioritize.hm(res, effect.size, qval.threshold)

    ret <- new("perturbation.hm.analysed",
           .gene.hits             =
             data.table::as.data.table(priorit$gene.hits),
           .gene.pathogen.hits    =
             data.table::as.data.table(priorit$gene.pathogen.hits),
           .gene.effects          = data.table::
             as.data.table(res$gene.effects),
           .gene.pathogen.effects =
             data.table::as.data.table(res$gene.pathogen.effects),
           .screen.type.effects   = data.table::
             as.data.table(res$infection.effects),
           .data                  = data.table::
             as.data.table(res$model.data@.data),
           .model.fit             = res$model,
           .is.bootstrapped       = res$btst,
           .params                = list(effect.size=effect.size,
                                         qval.threshold=qval.threshold))

    ret
  }
)

#' @noRd
#' @import data.table
#' @importFrom dplyr mutate full_join select group_by
.hm.model.data <- function(md, bootstrap.cnt)
{
  if (is.numeric(bootstrap.cnt) & bootstrap.cnt < 10 & bootstrap.cnt >= 1)
  {
    stop("Please use at least 10 bootstrap runs (better 100/1000).")
  }

  # save gene control mappings
  gene.control.map <-
    dplyr::select(md@.data, GeneSymbol, Control) %>%
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

  message("Fitting hm")
  fit.hm <- .hm(md)
  ref     <- .ranef(fit.hm)

  #calc fdrs
  gp.fdrs <- gp.fdrs(ref$gene.pathogen.effects)
  ge.fdrs <- ge.fdrs(md, ref, bootstrap.cnt)

  # set together the gene/fdr/effects and the mappings
  ges     <- dplyr::full_join(ref$gene.effects, gene.control.map,
                              by="GeneSymbol") %>%
    dplyr::full_join(dplyr::select(ge.fdrs$ret, GeneSymbol, Qval),
                     by="GeneSymbol")
  gps     <- dplyr::full_join(gp.fdrs$gene.pathogen.matrix,
                              gene.control.map, by="GeneSymbol")

  ret <- list(gene.effects=ges,
              gene.pathogen.effects=gps,
              screen.type.effects=ref$screen.type.effects,
              model.data=md,
              model=list(fit=fit.hm, gp.fdrs=gp.fdrs, ge.fdrs=ge.fdrs),
              btst=ge.fdrs$btst)
  ret
}

#' @noRd
.init.formula <- function()
{
  frm.str <- paste0("Readout ~ Condition + ",
                    "(1 | GeneSymbol) + (1 | Condition:GeneSymbol) + ",
                    "(1 | ScreenType) + (1 | Condition:ScreenType)")
  frm.str
}

#' @noRd
#' @import data.table
#' @importFrom lme4 lmer
#' @importFrom stats as.formula
.hm <- function(md)
{
  lme4::lmer(stats::as.formula(.init.formula()),
             data = md@.data, weights = md@.data$Weight, verbose = FALSE)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr mutate select
#' @importFrom lme4 ranef
.ranef <-  function(fit.hm)
{
  random.effects <- lme4::ranef(fit.hm)
  # create the data table with gene effects
  ge <- data.table::data.table(
    Effect     = random.effects[["GeneSymbol"]][,1],
    GeneSymbol = as.character(rownames(random.effects[["GeneSymbol"]]))
  )

  # create the data.table with gene-pathogen effects
  gpe <- data.table::data.table(
    gpe = random.effects[["Condition:GeneSymbol"]][,1],
    GenePathID =
      as.character(rownames(random.effects[["Condition:GeneSymbol"]]))) %>%
    dplyr::mutate(GeneSymbol = sub("^.+:", "", GenePathID))

  ie <- data.table::data.table(
    Effect = random.effects[["ScreenType"]][,1],
    ScreenType = as.character(rownames(random.effects[["ScreenType"]])))

  ga <- merge(gpe, ge, by = "GeneSymbol") %>%
    dplyr::mutate(Condition = sub(":.+$", "", GenePathID),
                  GeneConditionEffect = Effect + gpe) %>%
    dplyr::select(-GenePathID, -gpe, -Effect)

  list(gene.effects=ge,
       gene.pathogen.effects=ga,
       screen.type.effects=ie)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr filter group_by mutate
.prioritize.hm <- function(obj, eft, fdrt)
{
  ge <- dplyr::filter(obj$gene.effects, abs(Effect) >= eft)  %>%
    .[order(-abs(Effect))]
  if (obj$btst != 0)
    ge <- dplyr::filter(ge, Qval <= fdrt)
  gpe <- obj$gene.pathogen.effects %>%
    dplyr::filter(abs(Effect) >= eft, Qval <= fdrt)  %>%
    .[order(-abs(Effect))]

  list(gene.hits=ge, gene.pathogen.hits=gpe)
}
