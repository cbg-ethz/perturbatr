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

#' Off-target correct and rank by hit magnitude
#'
#' @export
#' @import data.table
#'
#' @param obj  a summarized data.table
#' @param path  path (or file) to the target-relation matrices
#' @param drop  drops all genes that are not found in every screen, i.e. a
#'  gene has to be found over all pathogens to be used
#' @param ... additional arguments
correct <- function(obj, path, drop, ...)
{
  UseMethod("correct")
}

#' @noRd
#' @export
#' @import data.table
correct.svd.data <- function(obj, path, drop, ...)
{
  if (missing(path))
      stop(paste("Please provide the path to your target-relation matrices",
                 "(or one of the files, wel'll parse it for you)!"))
  res <-  .off.target.correct(obj, path, drop, ...)
  class(res) <- c("svd.analysed.offtc", "svd.analysed", class(res))
  invisible(res)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr select
#' @importFrom dplyr filter
.off.target.correct <- function(obj, path, drop, ...)
{
  nor <- nrow(obj)
  # only take elements that have an entrez ID and a sirna ID
  res <-
    dplyr::filter(obj, !is.na(Entrez), !is.na(siRNAIDs)) %>%
    ungroup
  # drop genes that are not found in all screens
  if (drop)
  {
    vir.cnt <- length(unique(res$Virus))
    res <- dplyr::group_by(res, GeneSymbol) %>%
      dplyr::mutate(drop=(length(unique(Virus)) != vir.cnt)) %>%
      ungroup %>%
      dplyr::filter(!drop) %>% dplyr::select(-drop)
    message("Dropped rows with genes not found in every virus-screen!")
  }
  if (nrow(res) < nor) message("Rows with is.na(Entrez) have been removed!")
  res <- .gespeR(res, path)
  invisible(res)
}

#' @noRd
#' @import data.table
#' @import gespeR glmnet parallel doParallel Matrix
.gespeR <- function(obj, path)
{
  # load target relation matrix
  rel.mat <- .load.rds(path)
  # init the matrices so that all entrez id that are not in both matrices are discarded
  # and that both matrices have the right dimensional
  gesp.dat <- .init.frames(pheno.mat=obj, rel.mat=rel.mat)
  pheno.frame <- gesp.dat$pheno.frame
  rel.mat <- gesp.dat$rel.mat
  rm(gesp.dat)
  phenos <- gespeR::Phenotypes(Matrix::Matrix(pheno.frame$Readout, ncol=1),
                               ids=as.character(pheno.frame$siRNAIDs),
                               pnames="pheno",
                               type="SSP")
  # fit a model using gesper
  gesper.fit <- gespeR::gespeR(phenotypes=phenos,
                               target.relations=rel.mat,
                               mode = "cv",
                               ncores=1)
  gsps  <- as.vector(gesper.fit@GSP@values)
  gsp.idx <- which(!is.na(gsps))
  genes.g <- as.vector(gesper.fit@GSP@ids[gsp.idx])
  gps.g <- gsps[gsp.idx]

  fr <- dplyr::left_join(data.table(Entrez=as.integer(genes.g),Effect=gps.g),
                         dplyr::select(obj, GeneSymbol, Entrez) %>% unique,
                         by="Entrez") %>%
    dplyr::filter(!is.na(GeneSymbol))
  invisible(fr)
}

#' @import data.table
#' @import dtplyr
#' @importFrom dplyr select
#' @importFrom assertthat assert_that
.init.frames <- function(pheno.mat, rel.mat)
{
  # unnest dharmacon entries
  pheno.mat <- .preprocess.dharmacon.entries(pheno.mat) %>%
    dplyr::select(Virus, GeneSymbol, Entrez, siRNAIDs, Readout) %>%
    .[order(siRNAIDs)]
  # get names of sirnas and entrez ids
  pheno.readout  <- pheno.mat$Readout
  pheno.sirnas   <- as.character(pheno.mat$siRNAIDs)
  pheno.entrez   <- as.integer(pheno.mat$Entrez)
  rel.mat.sirnas <- as.character(rel.mat@siRNAs)
  rel.mat.entrez <- as.integer(rel.mat@genes)
  # check if matrix has correct names
  assertthat::assert_that(all(rel.mat.entrez == colnames(rel.mat@values)))
  assertthat::assert_that(all(rel.mat.sirnas == rownames(rel.mat@values)))
  # intercept of pheno and relmat sirna
  intr.sirnas        <- unique(intersect(pheno.sirnas, rel.mat.sirnas))
  # indexes if the relation matrix sirnas that we want to use
  rel.sirna.idxs     <- which(rel.mat.sirnas %in% intr.sirnas)
  pheno.sirna.idxs   <- which(pheno.sirnas %in% intr.sirnas)
  # remove indexes that are not found
  pheno.sirnas   <- pheno.sirnas[pheno.sirna.idxs]
  rel.mat.sirnas <- rel.mat.sirnas[rel.sirna.idxs]
  vals           <- rel.mat@values[rel.sirna.idxs, ]
  vals           <- vals[order(rownames(vals)), ]
  rel.mat.sirnas <- rownames(vals)
  # order the pheno sirnas according to the ordering of the rel.matrix sirnas
  # what and why:
  # 1) first check that the rel.mat.sirnas are UNIQUE
  # 2) match the pheno.sirnas to the rel.mat sirnas
  # (i.e. which index has a pheno sirnas in the rel.mat array,
  # these are M:1 mappings)
  # 3) get the indexes (order is need)
  assertthat::assert_that(length(rel.mat.sirnas) ==
                            length(unique(rel.mat.sirnas)))
  assertthat::assert_that(!is.unsorted(row.names(vals)))
  assertthat::assert_that(!is.unsorted(rel.mat.sirnas))
  assertthat::assert_that(!is.unsorted(pheno.sirnas))
  ord <- order(match(pheno.sirnas, rel.mat.sirnas))
  # rearrane the pheno.mat sirnas
  pheno.ordered <- pheno.sirnas[ord]
  #count how often th sirnas are availabel
  # ERROR: TABLE SORTS ALPHABETICALLY: REARRANGE HERE!!!!!!
  pheno.count   <- table(pheno.ordered)
  assertthat::assert_that(!is.unsorted(names(pheno.count)))
  assertthat::assert_that(!is.unsorted(names(pheno.count)))
  assertthat::assert_that(all(names(pheno.count) == rownames(vals)))
  # add redundant siRNAs to the target relation matrix
  vals <- vals[rep(1:nrow(vals), pheno.count), ]
  # assign to matrix
  rel.mat@siRNAs <- rownames(vals)
  rel.mat@values <- vals
  # create phenotype frame
  pheno.frame <- data.table(Entrez= pheno.mat$Entrez,
                            siRNAIDs=pheno.mat$siRNAIDs,
                            Readout=pheno.mat$Readout) %>%
    # take only sirnas that have been found
    .[pheno.sirna.idxs] %>%
    # reorder
    .[ord]
  #check if dimensions are so are right
  assertthat::assert_that(all(pheno.frame$siRNAIDs == rel.mat@siRNAs))
  assertthat::assert_that(all(pheno.frame$siRNAIDs ==
                                rownames(rel.mat@values)))
  # TODO some asserts if values are correct
  # THE VALUES ARE FALSE WHAT THE HELL
  list(pheno.frame=pheno.frame, rel.mat=rel.mat)
}

#' @import data.table
#' @importFrom dplyr select filter
#' @importFrom tidyr unnest
.preprocess.dharmacon.entries <- function(obj)
{
  obj.unnest <-
    obj %>%
    dplyr::mutate(siRNAIDs=strsplit(as.character(siRNAIDs), ",")) %>%
    tidyr::unnest(siRNAIDs)
  obj.unnest
}


#' @noRd
#' @importFrom gespeR TargetRelations
.load.rds <- function(path)
{
  rel.mat <- gespeR::TargetRelations(path)
  invisible(rel.mat)
}
