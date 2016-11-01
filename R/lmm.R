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
#' @param ...  additional parameters
#' \itemize{
#'  \item{ignore }{ remove sirnas that have been found less than \code{ignore} times}
#' }
lmm <- function(obj, drop=T, weights=NULL, rel.mat.path=NULL, ...)
  UseMethod("lmm", obj)

#' @noRd
#' @export
#' @import data.table
lmm.svd.data <- function(obj, drop=T, weights=NULL, rel.mat.path=NULL, ...)
{
  res  <- .lmm(obj, drop, weights, rel.mat.path)
  class(res) <- c("svd.analysed.pmm","svd.analysed", class(res))
  invisible(res)
}

#' @noRd
#' @import data.table
#' @import lme4
#' @importFrom dplyr mutate
#' @importFrom dplyr select
.lmm <- function(obj, drop, weights=NULL, rel.mat.path=NULL, ...)
{
  params <- base::list(...)
  ignore <- base::ifelse(hasArg(ignore) &&
                           is.numeric(params$ignore), params$ignore, 1)
  # init the data table for the LMM
  model.data <- .set.lmm.matrix(obj, drop, ignore, weights, rel.mat.path)
  # save gene control mappings
  gene.control.map <-
    dplyr::select(model.data, GeneSymbol, Control) %>%
    unique %>%
    dplyr::mutate(GeneSymbol=as.character(GeneSymbol))
  # fit the LMM
  fit.lmm <-
    lme4::lmer(Readout ~ Virus + (1 | GeneSymbol) + (1 | Virus:GeneSymbol),
               data = model.data, weights = model.data$Weight,
               verbose = FALSE)
  random.effects <- lme4::ranef(fit.lmm)
  # create the data table with gene effects
  ag <- data.table::data.table(
    Effect = random.effects[["GeneSymbol"]][,1],
    GeneSymbol = as.character(rownames(random.effects[["GeneSymbol"]])))
  # create the data.table with pathogen effects
  bcg <- data.table::data.table(
    bcg = random.effects[["Virus:GeneSymbol"]][,1],
    GenePathID = as.character(rownames(random.effects[["Virus:GeneSymbol"]]))) %>%
    dplyr::mutate(GeneSymbol = sub("^.+:", "", GenePathID))
  # create the table with gene-pathogen effects
  ccg <- base::merge(bcg, ag, by = "GeneSymbol") %>%
    dplyr::mutate(Virus = sub(":.+$", "", GenePathID),
                  GeneVirusEffect = Effect + bcg) %>%
    dplyr::select(-GenePathID, -bcg, -Effect)
  # calculate fdrs
  fdrs <- .fdr(ccg)
  # finalize output and return as list
  gene.effects <- dplyr::full_join(ag, gene.control.map,  by="GeneSymbol")
  gene.pathogen.effects <- dplyr::full_join(
    fdrs$gene.pathogen.matrix,
    gene.control.map,
    by="GeneSymbol"
  )
  ret <- base::list(
    gene.effects=gene.effects,
    gene.pathogen.effects=gene.pathogen.effects,
    model.data=model.data,
    fit=list(model=fit.lmm, fdrs=fdrs$fdrs))
  invisible(ret)
}

#' Create LMM matrix
#'
#' @noRd
#' @import data.table
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
.set.lmm.matrix <- function(obj, drop, ignore, weights=NULL, rel.mat.path=NULL)
{
  # set weights for the sirnas
  wei.dhar <- wei.am <- 1
  # check if weights are given
  if (!is.null(weights))
  {
    if (!is.null(weights$dharmacon))
    {
      wei.dhar <- weights$dharmacon
      cat(paste("Setting weights for Dharmacon library to", wei.dhar,"\n"))
    }
    if (!is.null(weights$ambion))
    {
      wei.am <- weights$ambion
      cat(paste("Setting weights for Ambion library to", wei.am,"\n"))
    }
  }
  # check if sirna-gene affinities are given
  else if (!is.null(rel.mat.path))
  {
    # TODO: this is the
    cat(paste("Setting weights for Dharmacon library to", wei.dhar,"\n"))
    # this i still have to implement
    rel.mat <- .load.rds(path)
  }
  # setup pmm data
  pmm.mat <-
    # subset the columns
    dplyr::select(obj, Entrez, GeneSymbol, Virus, Readout, Control, Library) %>%
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
  pmm.mat$Entrez     <- as.integer(pmm.mat$Entrez)
  pmm.mat$Virus      <- as.factor( pmm.mat$Virus)
  pmm.mat$VG         <- as.factor( pmm.mat$VG)
  pmm.mat$GeneSymbol <- as.factor( pmm.mat$GeneSymbol)
  pmm.mat$Weight     <- as.integer(pmm.mat$Weight)
  pmm.mat$Control    <- as.integer(pmm.mat$Control)
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
  invisible(pmm.mat)
}

#' @noRd
#' @importFrom stats reshape
#' @importFrom stats na.omit
.fdr <- function(obj, ...)
{
  # reshape the data matrix
  ccg.matrix <- stats::reshape(as.data.frame(obj), direction = "wide",
                               timevar = "Virus", idvar = "GeneSymbol")
  fdrs <- list()
  # calculate FDRs for every pathogen-gene effect
  for (i in unique(obj$Virus))
  {
    cl <- paste("GeneVirusEffect", i, sep = ".")
    fr <- paste("fdr", i, sep = ".")
    sub.ccg.matrix <- stats::na.omit(ccg.matrix[, c("GeneSymbol",cl) ])
    locf <- .localfdr(sub.ccg.matrix[, cl])
    fdrs[[i]] <- locf
    sub.ccg.matrix <- cbind(sub.ccg.matrix, fdr = locf$fdr)
    nmiss <- sum(is.na(ccg.matrix[, cl]))
    dmiss <- ccg.matrix[is.na(ccg.matrix[, cl ]), c("GeneSymbol", cl)]
    sub.ccg.matrix <- rbind(sub.ccg.matrix, cbind(dmiss, fdr = rep(NA, nmiss)))
    colnames(sub.ccg.matrix)[colnames(sub.ccg.matrix) == "fdr"] <- fr
    ccg.matrix <- merge(ccg.matrix, sub.ccg.matrix[, c("GeneSymbol", fr)],
                        by = "GeneSymbol")
  }
  # parse the result matrix into a gene matrix
  # TODO: extra function
  gene.effect.mat <-
    ccg.matrix[,c(1, grep("GeneVirusEffect", colnames(ccg.matrix)))]
  colnames(gene.effect.mat) <-
    sub("GeneVirusEffect.", "", colnames(gene.effect.mat))
  gene.effect.mat <- gene.effect.mat %>%
    tidyr::gather(Virus, Effect, 2:ncol(gene.effect.mat)) %>%
    as.data.table
  # parse the result matrix into a fdr matrix
  fdr.mat         <- ccg.matrix[,c(1, grep("fdr", colnames(ccg.matrix)))]
  colnames(fdr.mat) <- sub("fdr.", "", colnames(fdr.mat))
  fdr.mat <- fdr.mat %>%
    tidyr::gather(Virus, FDR, 2:ncol(fdr.mat)) %>%
    as.data.table
  # join both matrices
  gene.pathogen.matrix <-
    dplyr::full_join(gene.effect.mat, fdr.mat, by=c("GeneSymbol", "Virus"))
  list(ccg.matrix=ccg.matrix,
       fdrs=fdrs,
       gene.pathogen.matrix=gene.pathogen.matrix )
}
