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
#' @param loocv  the number of loocv runs you want to do in order to estimate a significance level for the gene effects
#' @param ...  additional parameters
#' \itemize{
#'  \item{ignore }{ remove sirnas that have been found less than \code{ignore} times}
#' }
lmm <- function(obj, drop=T,
                weights=NULL, rel.mat.path=NULL,
                loocv=5, ...)
{
  UseMethod("lmm")
}

#' @export
#' @import data.table
#' @method lmm svd.data
lmm.svd.data <- function(obj, drop=T,
                         weights=NULL, rel.mat.path=NULL,
                         loocv=5, ...)
{
  res  <- .lmm.svd.data(obj, drop,
                        weights, rel.mat.path, loocv, ...)
  class(res) <- c("svd.analysed.pmm","svd.analysed", class(res))
  invisible(res)
}

#' @export
#' @import data.table
#' @method lmm svd.lmm.model.data
lmm.svd.lmm.model.data <- function(obj, drop=T,
                                   weights=NULL, rel.mat.path=NULL,
                                   loocv, ...)
{
  res <- .lmm.model.data(obj, loocv)
  class(res) <- c("svd.analysed.pmm","svd.analysed", class(res))
  invisible(res)
}

#' @noRd
#' @import data.table
#' @importFrom methods hasArg
.lmm.svd.data <- function(obj, drop,
                          weights=NULL, rel.mat.path=NULL,
                          loocv, ...)
{
  params <- list(...)
  ignore <- ifelse(methods::hasArg(ignore) && is.numeric(params$ignore),
                   params$ignore, 1)
  # init the data table for the LMM
  md <- .set.lmm.matrix(obj, drop, ignore, weights, rel.mat.path)
  .lmm.model.data(md, loocv)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr mutate full_join select group_by
.lmm.model.data <- function(md, loocv)
{
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
  # finalize output and return as list
  message("LOOCV for significance estimation")
  ge.fdrs <- .lmm.significant.hits(md)
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
#' @import data.table
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
.set.lmm.matrix <- function(obj, drop, ignore, weights=NULL, rel.mat.path=NULL)
{
  # TODO: change this to pooled and unpooled
  # set weights for the sirnas
  wei.dhar <- wei.am <- 1
  # check if weights are given
  if (!is.null(weights))
  {
    if (!is.null(weights$dharmacon))
    {
      wei.dhar <- weights$dharmacon
      message(paste("Setting weights for Dharmacon library to", wei.dhar,"\n"))
    }
    if (!is.null(weights$ambion))
    {
      wei.am <- weights$ambion
      message(paste("Setting weights for Ambion library to", wei.am,"\n"))
    }
  }
  # check if sirna-gene affinities are given
  else if (!is.null(rel.mat.path))
  {
    # TODO: this is the
    cat(paste("Setting weights for Dharmacon library to", wei.dhar,"\n"))
    # this i still have to implement
    stop("to do")
    rel.mat <- .load.rds(path)
  }
  # setup pmm data
  pmm.mat <-
    # subset the columns
    dplyr::select(obj, Entrez, GeneSymbol, Virus, Readout, Control, Library,
                  Cell, InfectionType, Design, ReadoutType) %>%
    # only take entries that have a genesymbol
    dplyr::filter(!is.na(GeneSymbol)) %>%
    # dont take positive controls since these are different between the pathonges
    # negative control on the other hand should be fine
    dplyr::filter(Control != 1)
  data.table::setDT(pmm.mat)[Library == "Dharmacon", Weight := wei.dhar]
  data.table::setDT(pmm.mat)[Library == "Ambion",    Weight := wei.am]
  # set a column that concats virus and genesymbol
  data.table::setDT(pmm.mat)[, VG := paste(Virus, GeneSymbol, sep=":")]
  #  remove librarz
  data.table::setDT(pmm.mat)[, Library := NULL]
  ## cast some columns
  pmm.mat$Entrez        <- as.integer(pmm.mat$Entrez)
  pmm.mat$Virus         <- as.factor(pmm.mat$Virus)
  pmm.mat$VG            <- as.factor(pmm.mat$VG)
  pmm.mat$Cell          <- as.factor(pmm.mat$Cell)
  pmm.mat$ReadoutType   <- as.factor(pmm.mat$ReadoutType)
  pmm.mat$InfectionType <- as.factor(pmm.mat$InfectionType)
  pmm.mat$Design        <- as.factor(pmm.mat$Design)
  pmm.mat$GeneSymbol    <- as.factor(pmm.mat$GeneSymbol)
  pmm.mat$Weight        <- as.double(pmm.mat$Weight)
  pmm.mat$Control       <- as.integer(pmm.mat$Control)
  pmm.mat <-
    # remove entries that have nan as readout
    dplyr::filter(pmm.mat, !is.na(Readout)) %>%
    dplyr::group_by(VG) %>%
    # count how often VG is in the data
    dplyr::mutate(cnt=n()) %>%
    ungroup %>%
    # throw away VGs that are less than ignore
    dplyr::filter(cnt >= ignore) %>%
    # remove cnt column
    dplyr::select(-cnt)
  # drop genes that are not found in ALL pathogens
  if (drop)
  {
    # count how many viruses are in the date sets
    vir.cnt <- length(unique(pmm.mat$Virus))
    pmm.mat <-
      # group by genesymbol
      dplyr::group_by(pmm.mat, GeneSymbol) %>%
      # count if the genes are in all viruses
      # and compare if it matches the virus count
      dplyr::mutate(drop=(length(unique(Virus)) != vir.cnt)) %>%
      ungroup %>%
      # remove genes that are not in all genes
      dplyr::filter(!drop) %>%
      # remove drop column
      dplyr::select(-drop)
  }
  pmm.mat <- droplevels(pmm.mat)
  class(pmm.mat) <- c("svd.lmm.model.data", class(pmm.mat))
  pmm.mat
}

#' @noRd
.init.formula <- function()
{
  frm.str <- paste0("Readout ~ Virus + ",
                    "(1 | GeneSymbol) + (1 | Virus:GeneSymbol) + ",
                    "(1 | InfectionType) + (1 | Virus:InfectionType)")
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
.lmm.significant.hits <- function(model.data, padj=c("bonf", "BH"))
{
  padj <- match.arg(padj)
  # do 5 loocv iterations
  li <- list()
  i <- 1
  ctr <- 1
  repeat
  {
    ctr <- ctr + 1
    tryCatch({
      loocv.sample <- as.svd.lmm.model.data(loocv(model.data, i))
      # here use the lmer params
      lmm.fit <- .lmm(loocv.sample)
      re <- .ranef(lmm.fit)
      da <- data.table::data.table(loocv=paste0("LOOCV_", i),
                                   Effect=re$gene.effects$Effect,
                                   GeneSymbol=re$gene.effects$GeneSymbol)
      li[[i]] <- da
      i <- i + 1
    }, error=function(e) { print(paste("Didn't fit:", i, ", error:", e)); i <<- 1000 },
    warning=function(e) { print(e)})
    if (i > 5)
      break
    if (ctr == 100)
      stop("Breaking after 100 mis-trials!")
  }
  flat.dat <- do.call("rbind", lapply(li, function(e) e))
  dat <-  flat.dat %>%
    dplyr::group_by(GeneSymbol) %>%
    dplyr::summarise(MeanLOOCVEffect=mean(Effect, na.rm=T),
                     Pval=.ttest(GeneSymbol, Effect, 0)) %>%
    ungroup %>%
    dplyr::mutate(FDR=p.adjust(Pval, method=padj)) %>%
    .[order(FDR)]
  dat <- dplyr::left_join(dat, tidyr::spread(flat.dat, loocv, Effect), bby="GeneSymbol")
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
    Effect = random.effects[["InfectionType"]][,1],
    InfectionType = as.character(rownames(random.effects[["InfectionType"]])))
  # table for virus-infection types
  # ipe <- data.table::data.table(
  #   ie = random.effects[["Virus:InfectionType"]][,1],
  #   InfectionPathID = as.character(rownames(random.effects[["Virus:InfectionType"]])))
  # create the table with gene-pathogen effects
  ga <- base::merge(gpe, ge, by = "GeneSymbol") %>%
    dplyr::mutate(Virus = sub(":.+$", "", GenePathID),
                  GeneVirusEffect = Effect + gpe) %>%
    dplyr::select(-GenePathID, -gpe, -Effect)
  list(gene.effects=ge,
       gene.pathogen.effects=ga,
       infection.effects=ie)
}
