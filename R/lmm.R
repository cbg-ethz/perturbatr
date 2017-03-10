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


#' Fit an LMM to the data and calculate local false discovery rates.
#'
#' @export
#'
#' @import data.table
#'
#' @param obj  an svd.data object
#' @param drop  boolean flag if all entries should be dropped that are not found in every virus
#' @param weights a list of weights
#' @param rel.mat.path  the (optional) path to a target relation matrix that is going to be used for
#' @param bootstrap.cnt  the number of loocv runs you want to do in order to estimate a significance level for the gene effects
#' @param ignore  ignore siRNAs that have been seen only once per group
#' @param ...  additional parameters
lmm <- function(obj, drop=T,
                weights=NULL, rel.mat.path=NULL,
                bootstrap.cnt=F, ignore=1, ...)
{
  UseMethod("lmm")
}

#' @export
#' @import data.table
#' @method lmm svd.data
lmm.svd.data <- function(obj, drop=T,
                         weights=NULL, rel.mat.path=NULL,
                         bootstrap.cnt=F, ignore=1, ...)
{
  res  <- .lmm.svd.data(obj, drop,
                        weights, rel.mat.path, bootstrap.cnt, ignore, ...)
  class(res) <- c("svd.analysed.pmm","svd.analysed", class(res))
  invisible(res)
}

#' @export
#' @import data.table
#' @method lmm svd.lmm.model.data
lmm.svd.lmm.model.data <- function(obj, drop=T,
                                   weights=NULL, rel.mat.path=NULL,
                                   bootstrap.cnt=F, ignore=1, ...)
{
  res <- .lmm.model.data(obj, bootstrap.cnt)
  class(res) <- c("svd.analysed.pmm","svd.analysed", class(res))
  invisible(res)
}

#' @noRd
#' @import data.table
#' @importFrom methods hasArg
.lmm.svd.data <- function(obj, drop,
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
#' @importFrom assertthat assert_that
.weights <- function(obj, weights, rel.mat.path)
{
  if (is.null(weights)) return(1)
  else if (!is.list(weights)) stop("Please give a list argument")
  els <- names(weights)
  ret <- rep(1, nrow(obj))
  for (el in els)
  {
    idxs <- switch(el,
                   "pooled"  =which(obj$Design == "pooled"),
                   "single"=which(obj$Design == "single"),
                   stop("Please provide 'single'/'pooled' list names for setting weights"))
    ret[idxs] <- as.numeric(weights[[el]])
    message(paste("Setting", length(idxs), el,  "well weights to:", as.numeric(weights[[el]])))
  }
  assertthat::assert_that(length(ret) == nrow(obj))
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
#' @importFrom dplyr mutate select left_join
#' @importFrom tidyr spread
#' @importFrom lme4 ranef
.lmm.significant.hits <- function(model.data, bootstrap.cnt, padj=c("BH", "bonf"))
{
  padj <- match.arg(padj)
  li <- list()
  i <- 1
  ctr <- 1
  mistrial.cnt <- bootstrap.cnt * 10
  repeat
  {
    ctr <- ctr + 1
    tryCatch({
      bt.sample <- as.svd.lmm.model.data(bootstrap(model.data))
      # here use the lmer params
      lmm.fit <- .lmm(bt.sample)
      re <- .ranef(lmm.fit)
      da <- data.table::data.table(bootstrap=paste0("Bootstrap_", sprintf("%03i", i)),
                                   Effect=re$gene.effects$Effect,
                                   GeneSymbol=re$gene.effects$GeneSymbol)
      li[[i]] <- da
      i <- i + 1
    }, error=function(e) { print(paste("Didn't fit:", i, ", error:", e)); i <<- 1000 },
    warning=function(e) { print(e)})
    if (i > bootstrap.cnt)
      break
    if (ctr == mistrial.cnt)
      stop(paste0("Breaking after ", mistrial.cnt ," mis-trials!"))
  }
  flat.dat <- do.call("rbind", lapply(li, function(e) e))
  # TODO modify that accordingly
  dat <-  flat.dat %>%
    dplyr::group_by(GeneSymbol) %>%
    dplyr::summarise(MeanBootstrap=mean(Effect, na.rm=T),
                     Pval=.ttest(GeneSymbol, Effect, 0)) %>%
    ungroup %>%
    dplyr::mutate(FDR=p.adjust(Pval, method=padj)) %>%
    .[order(FDR)]
  dat <- dplyr::left_join(dat, tidyr::spread(flat.dat, bootstrap, Effect), by="GeneSymbol")
  dat
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
