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
#' @include inference_lmm_fdr.R
#' @include inference_lmm_locfdr.R


#' @title Jointly analyse multiple genetic perturbation screens using a
#'  hierchical model
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
#' @param obj  an \code{PerturbationData} object
#' @param formula  a \code{formula} object that is used to model the readout
#'  of your data set. If no formula is provided, the formula
#'  `Readout ~ Condition + (1|GeneSymbol) +  (1|Condition:GeneSymbol)` is used. For
#'  other data sets with more variables, it might makes sense to use other
#'  fixed and random effects
#' @param drop  boolean if genes that are not found in every Condition should
#'  be dropped
#' @param weights  a numeric vector used as weights for the single
#'  perturbations
#' @param bootstrap.cnt  the number of bootstrap runs you want to do in order
#'  to estimate a significance level for the gene effects
#' @param effect.size  the absolute, relative effect size used for
#'  hit prioritization
#' @param qval.threshold  the q-value threshold used for hit prioritization
#'  if \code{bootstrap.cnt} is set
#'
#' @return returns a \code{HMAnalysedPerturbationData} object
#'
#' @examples
#'  data(rnaiscreen)
#'  rnaiscreen.normalized <- preprocess(rnaiscreen, normalize="robust-z.score")
#'  res                   <- hm(rnaiscreen.normalized, effect.size=0.01)
setGeneric(
  "hm",
  function(obj,
           formula=Readout ~ Condition+(1|GeneSymbol)+(1|Condition:GeneSymbol),
           drop=TRUE,
           weights=1,
           bootstrap.cnt=0,
           effect.size=0.05,
           qval.threshold=.2)
  {
    standardGeneric("hm")
  }
)


#' @rdname hm-methods
#' @aliases hm,PerturbationData-method
#' @import data.table
setMethod(
  "hm",
  signature = signature(obj="PerturbationData"),
  function(obj,
           formula=Readout ~ Condition+(1|GeneSymbol)+(1|Condition:GeneSymbol),
           drop=TRUE,
           weights=1,
           bootstrap.cnt=0,
           effect.size=0.05,
           qval.threshold=.2)
  {
    md <- setModelData(obj, drop, weights)
    .hm(md, formula, drop, weights, bootstrap.cnt, effect.size, qval.threshold)
  }
)


#' @noRd
#' @import data.table
#' @importFrom methods new
.hm  <- function(
    obj,
    formula=Readout ~ Condition+(1|GeneSymbol)+(1|Condition:GeneSymbol),
    drop=TRUE,
    weights=1,
    bootstrap.cnt=0,
    effect.size=0.01,
    qval.threshold=.2)
{
  res     <- .hm.model.data(obj, bootstrap.cnt)
  priorit <- .prioritize.hm(res, effect.size, qval.threshold)

  ret <- methods::new("HMAnalysedPerturbationData",
          geneHits = data.table::as.data.table(priorit$gene.hits),
          nestedGeneHits = data.table::as.data.table(priorit$nested.gene.hits),
          geneEffects = data.table::as.data.table(res$gene.effects),
          nestedGeneEffects = data.table::as.data.table(res$nested.gene.effects),
          data = data.table::as.data.table(res$model.data@.data),
          modelFit = res$model,
          isBootstrapped = res$btst,
          params = list(formula=formula,
                        effect.size=effect.size,
                        qval.threshold=qval.threshold,
                        weights=weights,
                        drop=drop,
                        bootstrap.cnt=bootstrap.cnt))

  ret
}


#' @noRd
#' @import data.table
#' @importFrom dplyr mutate full_join select group_by n pull
.hm.model.data <- function(md, bootstrap.cnt)
{
  if (is.numeric(bootstrap.cnt) & bootstrap.cnt < 10 & bootstrap.cnt >= 1)
  {
    stop("Please use at least 10 bootstrap runs (better 100/1000).")
  }

  # save gene control mappings
  gene.control.map <-
    dplyr::select(md, GeneSymbol, Control) %>%
    unique() %>%
    dplyr::mutate(GeneSymbol=as.character(GeneSymbol))
  mult.gen.cnt <- gene.control.map %>%
                     dplyr::group_by(GeneSymbol) %>%
                     dplyr::mutate(cnt = n()) %>%
    dplur::pull(cnt) %>%
    unique()
  if (length(mult.gen.cnt) != 1)
  {
    warning(paste("Found multiple gene-control entries for several genes,",
                   "i.e. several genes are both control and not control!"))
  }

  message("Fitting hierarchical model")
  fit.hm   <- .hm(md, formula)
  ref      <- .ranef(fit.hm)
  ge.fdrs.ef  <- ge.fdrs(md, ref, bootstrap.cnt)
  nge.fdrs.ef <- nge.fdrs(ref$nested.gene.effects)

  # set together the gene/fdr/effects and the mappings
  ges <- dplyr::full_join(ref$gene.effects, gene.control.map,
                          by="GeneSymbol") %>%
    dplyr::full_join(dplyr::select(ge.fdrs.ef$ret, GeneSymbol, Qval),
                                   by="GeneSymbol")
  nges <- dplyr::full_join(nge.fdrs.ef$nested.gene.matrix,
                           gene.control.map, by="GeneSymbol")

  ret <- list(gene.effects        = ges,
              nested.gene.effects = nges,
              model.data          = md,
              model = list(fit=fit.hm, ge.fdrs=ge.fdrs.ef, nge.fdrs=nge.fdrs.ef),
              btst                = ge.fdrs$btst)
  ret
}


#' @noRd
#' @import data.table
#' @importFrom lme4 lmer
#' @importFrom stats as.formula
.hm <- function(md, formula)
{
  lme4::lmer(stats::as.formula(formula),
             data = md, weights = md$Weight, verbose = FALSE)
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
    GeneSymbol = as.character(rownames(random.effects[["GeneSymbol"]])))

  # create the data.table with gene-pathogen effects
  gpe <- data.table::data.table(
    gpe = random.effects[["Condition:GeneSymbol"]][,1],
    GenePathID =
      as.character(rownames(random.effects[["Condition:GeneSymbol"]]))) %>%
    dplyr::mutate(GeneSymbol = sub("^.+:", "", GenePathID))

  ga <- merge(gpe, ge, by = "GeneSymbol") %>%
    dplyr::mutate(Condition = sub(":.+$", "", GenePathID),
                  GeneConditionEffect = Effect + gpe) %>%
    dplyr::select(-GenePathID, -gpe, -Effect)

  list(gene.effects        = ge,
       nested.gene.effects = ga)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr filter group_by mutate
.prioritize.hm <- function(obj, eft, fdrt)
{
  ge <- dplyr::filter(obj$gene.effects, abs(Effect) >= eft)  %>%
    .[order(-abs(Effect))]

  if (obj$btst)
    ge <- dplyr::filter(ge, Qval <= fdrt)

  nge <- obj$nested.gene.effects %>%
    dplyr::filter(abs(Effect) >= eft, Qval <= fdrt) %>%
    dplyr::arrange(-abs(Effect))


  list(gene.hits=ge, nested.gene.hits=nge)
}
