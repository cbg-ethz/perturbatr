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
#'  hierarchical model
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
#' @import tibble
#' @import formula.tools
#'
#' @param obj  an \code{PerturbationData} object
#' @param formula  a \code{formula} object that is used to model the readout
#'  of your data set. If no formula is provided, the formula
#'  `Readout ~ Condition + (1|GeneSymbol) +  (1|Condition:GeneSymbol)` is used.
#' For other data sets with more variables, it might makes sense to use other
#'  fixed and random effects
#' @param drop  boolean if genes that are not found in every Condition should
#'  be dropped
#' @param weights  a numeric vector used as weights for the single
#'  perturbations
#' @param bootstrap.cnt  the number of bootstrap runs you want to do in order
#'  to estimate a significance level for the gene effects
#'
#' @return returns a \code{HMAnalysedPerturbationData} object
#'
#' @examples
#'  data(rnaiscreen)
#'  res <- hm(rnaiscreen)
setGeneric(
  "hm",
  function(obj,
           formula=Readout ~ Condition+(1|GeneSymbol)+(1|Condition:GeneSymbol),
           drop=TRUE,
           weights=1,
           bootstrap.cnt=0)
  {
    standardGeneric("hm")
  }
)


#' @rdname hm-methods
#' @aliases hm,PerturbationData-method
#' @import tibble

setMethod(
  "hm",
  signature = signature(obj="PerturbationData"),
  function(obj,
           formula=Readout ~ Condition+(1|GeneSymbol)+(1|Condition:GeneSymbol),
           drop=TRUE,
           weights=1,
           bootstrap.cnt=0)
  {
    .check.formula(formula)
    md <- setModelData(obj, drop, weights)
    .hm(md, as.character(formula), drop, weights, bootstrap.cnt)
  }
)

#' @noRd
#' @importFrom formula.tools rhs.vars
.check.formula <- function(formula)
{

  needed.frms <- c("^Condition$",
                   "^1 \\| GeneSymbol$",
                   "^1 \\| GeneSymbol:Condition|Condition:GeneSymbol$")
  txt         <- c("fixed effect: 'Condition'",
                   "random intercept: '1 | GeneSymbol'",
                   "random intercept: '1 | Condition:GeneSymbol'")
  pre <- "Plase build your formula with a"
  frm.rhs <- formula.tools::rhs.vars(formula)
  for (i in seq(needed.frms)) {
    if (is.na(grep(needed.frms[i], frm.rhs)[1]))
      stop(paste(pre, txt[i]))
  }
}


#' @noRd
#' @import tibble
#' @importFrom methods new
.hm  <- function(obj, formula, drop, weights, bootstrap.cnt)
{
  res     <- .hm.model.data(obj, formula, bootstrap.cnt)
  params  <- list(formula = formula,
                  weights = weights,
                  drop    = drop,
                  bootstrap.cnt = bootstrap.cnt)

  ret <- methods::new("HMAnalysedPerturbationData",
          geneEffects = tibble::as.tibble(res$gene.effects),
          nestedGeneEffects = tibble::as.tibble(res$nested.gene.effects),
          dataSet = tibble::as.tibble(res$model.data),
          modelFit = res$model,
          isBootstrapped = res$btst,
          params = params)

  ret
}


#' @noRd
#' @import tibble
#' @importFrom dplyr mutate full_join select group_by n pull
#' @importFrom rlang .data
.hm.model.data <- function(md, formula, bootstrap.cnt)
{
  if (is.numeric(bootstrap.cnt) & bootstrap.cnt < 10 & bootstrap.cnt >= 1) {
    stop("Please use at least 10 bootstrap runs (better 100/1000).")
  }

  # save gene control mappings
  gene.control.map <- unique(
    dplyr::select(md, .data$GeneSymbol, .data$Control))
  gene.control.map <- dplyr::mutate(
    gene.control.map, "GeneSymbol" = as.character(.data$GeneSymbol))
  mult.gen.cnt <- dplyr::group_by(gene.control.map, .data$GeneSymbol)
  mult.gen.cnt <- dplyr::mutate(mult.gen.cnt, "cnt" = n())
  mult.gen.cnt <- unique(dplyr::pull(mult.gen.cnt, .data$cnt))

  if (length(mult.gen.cnt) != 1) {
    warning(paste("Found multiple gene-control entries for several genes,",
                   "i.e. several genes are both control and not control!"))
  }

  message("Fitting hierarchical model")
  fit.hm      <- .hm.fit(md, formula)
  ref         <- .ranef(fit.hm)
  ge.fdrs.ef  <- ge.fdrs(md, ref, bootstrap.cnt, formula)
  nge.fdrs.ef <- nge.fdrs(ref$nested.gene.effects)

  # set together the gene/fdr/effects and the mappings
  ges <- dplyr::full_join(ref$gene.effects, gene.control.map, by="GeneSymbol")
  ges <- dplyr::full_join(
    ges,
    dplyr::select(ge.fdrs.ef$ret, .data$GeneSymbol, .data$Qval),
    by="GeneSymbol"
  )
  nges <- dplyr::full_join(
    nge.fdrs.ef$nested.gene.matrix,
    gene.control.map,
    by="GeneSymbol"
  )

  ret <- list(gene.effects        = ges,
              nested.gene.effects = nges,
              model.data          = md,
              model = list(fit=fit.hm,ge.fdrs=ge.fdrs.ef, nge.fdrs=nge.fdrs.ef),
              btst                = ge.fdrs.ef$btst)
  ret
}


#' @noRd
#' @import tibble
#' @importFrom lme4 lmer
#' @importFrom stats as.formula
.hm.fit <- function(md, formula)
{
  lme4::lmer(stats::as.formula(formula),
             data = md, weights = md$Weight, verbose = FALSE)
}

#' @noRd
#' @import tibble
#' @importFrom dplyr mutate select
#' @importFrom lme4 ranef
#' @importFrom rlang .data
.ranef <-  function(fit.hm)
{
  random.effects <- lme4::ranef(fit.hm)
  # create the data table with gene effects
  ge <- tibble::tibble(
    Effect     = random.effects[["GeneSymbol"]][,1],
    GeneSymbol = as.character(rownames(random.effects[["GeneSymbol"]])))

  # create the tibble with gene-pathogen effects
  gpe <- tibble::tibble(
    gpe = random.effects[["Condition:GeneSymbol"]][,1],
    GenePathID =
      as.character(rownames(random.effects[["Condition:GeneSymbol"]]))) %>%
    dplyr::mutate("GeneSymbol" = sub("^.+:", "", .data$GenePathID))

  ga <- base::merge(gpe, ge, by = "GeneSymbol") %>%
    dplyr::mutate("Condition" = sub(":.+$", "", .data$GenePathID),
                  "GeneConditionEffect" = .data$Effect + .data$gpe) %>%
    dplyr::select(-.data$GenePathID, -.data$gpe, -.data$Effect)

  list(gene.effects        = ge,
       nested.gene.effects = ga)
}
