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

#' Do overrepresentation analysis based on the hypergeometric test
#'
#' @export
#' @docType methods
#' @rdname ora-methods
#'
#' @param hit.list  a list of identified hit genes
#' @param universe  a list of all genes
#' @param db  either run ORA on KEGG or on Gene Ontology
#' @param ...  additional parameters
setGeneric(
  "ora",
   function(hit.list, universe, db=c("kegg", "go") , ...) standardGeneric("ora")
)

#' @rdname ora-methods
#' @aliases ora,character,character-method
setMethod(
  "ora",
  c(hit.list="character", universe="character"),
  function(hit.list, universe, db=c("kegg", "go"), ...)
    .ora(hit.list, universe, db, ...)
)

#' @noRd
#' @importFrom GOstats hyperGTest
.ora <- function(hit.list, universe, db=c("kegg", "go"), ...)
{
  db <- match.arg(db)

  hit.list <- .to.entrez(hit.list)
  universe <- .to.entrez(universe)
  message(paste("Doing ORA on:", db))
  dat <- switch(db,
                "kegg"=.gsea.kegg(hit.list, universe, ...),
                "go"  =.gsea.go(hit.list, universe, ...))
  ora <- GOstats::hyperGTest(dat)
  # this is tedious. they should correct this
  test.count <- switch(db,
                       "go"=length(ora@pvalue.order),
                       "kegg"=length(ora@pvalues))
  summ <- summary(ora)
  summ$Qvalue <- summ$Pvalue * test.count
  li <- list(ora=ora, summary=summ)
  class(li) <- c(class(li), "svd.ora")
  invisible(li)
}

#' @noRd
#' @importFrom methods hasArg
.gsea.go <- function(hit.list, universe, ...)
{
  params <- list(...)
  ontology <- ifelse(methods::hasArg("ontology"), params$ontology, "BP")
  p.value <- ifelse(methods::hasArg("p.value"), params$p.value, .05)
  GOparams <- new("GOHyperGParams",
                  geneIds = unique(hit.list$gene_id),
                  universeGeneIds = unique(universe$gene_id),
                  annotation="hgu95av2.db",
                  ontology=ontology,
                  pvalueCutoff=p.value,
                  conditional=TRUE,
                  testDirection="over")
  invisible(GOparams)
}

#' @noRd
#' @importFrom methods hasArg
.gsea.kegg <- function(hit.list, universe, ...)
{
  params <- list(...)
  p.value <- ifelse(methods::hasArg("p.value"), params$p.value, .05)
  KEGGparams <- new("KEGGHyperGParams",
                  geneIds =  unique(hit.list$gene_id),
                  universeGeneIds =  unique(universe$gene_id),
                  annotation="hgu95av2.db",
                  pvalueCutoff=p.value,
                  testDirection="over")
  invisible(KEGGparams)
}

#' @noRd
go.mapping <- function(dat)
{
  UseMethod("go.mapping", dat)
}

#' @noRd
#'
#' @import data.table
#' @import org.Hs.eg.db
#'
#' @importFrom AnnotationDbi toTable
#' @importFrom dplyr filter
go.mapping.integer <- function(dat)
{
  frame <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egGO)
  go.frame.data <- data.table::data.table(GO=frame$go_id,
                                          Entrez=frame$gene_id)
  res <- unique(dplyr::filter(go.frame.data, Entrez %in% dat)$GO)
  res
}

#' @noRd
kegg.mapping <- function(dat)
{
    UseMethod("kegg.mapping", dat)
}

#' @noRd
#'
#' @import data.table
#' @import org.Hs.eg.db
#' @importFrom AnnotationDbi toTable
#' @importFrom dplyr filter left_join
.to.entrez <-function(dat)
{
  frame.hugo <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egSYMBOL) %>%
    as.data.table
  dat <- dplyr::left_join(data.table(symbol=toupper(dat)), frame.hugo, by="symbol")
  dat
}

#' @noRd
#'
#' @import data.table
#' @import org.Hs.eg.db
#' @importFrom AnnotationDbi toTable
#' @importFrom dplyr filter
kegg.mapping.integer <-
function(dat)
{
  frame <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egPATH)
  kegg.frame.data <- data.table::data.table(Kegg=frame$path_id,
                                            Entrez=frame$gene_id)
  res <- unique(dplyr::filter(kegg.frame.data, Entrez %in% dat)$Kegg)
  res
}
