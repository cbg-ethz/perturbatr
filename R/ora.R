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
#' @aliases ora,integer,integer-method
setMethod(
  "ora",
  c(hit.list="integer", universe="integer"),
  function(hit.list, universe, db=c("kegg", "go"), ...)
    ora.default(hit.list, universe, db, ...)
)

#' @rdname ora-methods
#' @aliases ora,character,character-method
setMethod(
  "ora",
  c(hit.list="character", universe="character"),
  function(hit.list, universe, db=c("kegg", "go"), ...)
    ora.default(hit.list, universe, db, ...)
)

#' @noRd
#' @importFrom GOstats hyperGTest
ora.default <-
function
(
  hit.list,
  universe,
  db=c("kegg", "go"),
  ...
)
{
  db <- match.arg(db)
  message(paste("Doing ORA on:", db))
  dat <- switch(db,
                "kegg"=.gsea.kegg(hit.list, universe, ...),
                "go"=.gsea.go(hit.list, universe, ...))
  ora <- GOstats::hyperGTest(dat)
  summ <- summary(ora)
  li <- list(ora=ora, summary=summ)
  class(li) <- c(class(li), "svd.ora")
  invisible(li)
}

#' @noRd
.gsea.go <-
function
(
  hit.list,
  universe,
  ...
)
{
  params <- list(...)
  ontology <- ifelse(hasArg("ontology"), params$ontology, "BP")
  p.value <- ifelse(hasArg("p.value"), params$p.value, .05)
  GOparams <- new("GOHyperGParams",
                  geneIds = hit.list,
                  universeGeneIds = universe,
                  annotation="hgu95av2.db",
                  ontology=ontology,
                  pvalueCutoff=p.value,
                  conditional=TRUE,
                  testDirection="over")
  invisible(GOparams)
}

#' @noRd
.gsea.kegg <-
  function
(
  hit.list,
  universe,
  ...
)
{
  params <- list(...)
  p.value <- ifelse(hasArg("p.value"), params$p.value, .05)
  KEGGparams <- new("KEGGHyperGParams",
                  geneIds = hit.list,
                  universeGeneIds = universe,
                  annotation="hgu95av2.db",
                  pvalueCutoff=p.value,
                  testDirection="over")
  invisible(KEGGparams)
}


go.mapping <-
function(dat)
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
go.mapping.integer <-
function(dat)
{
  frame <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egGO)
  go.frame.data <- data.table::data.table(GO=frame$go_id,
                                          Entrez=frame$gene_id)
  res <- unique(dplyr::filter(go.frame.data, Entrez %in% dat)$GO)
  res
}

#' @noRd
kegg.mapping <-
function(dat)
{
    UseMethod("kegg.mapping", dat)
}

#' @noRd
#'
#' @import data.table
#' @import org.Hs.eg.db
#'
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
