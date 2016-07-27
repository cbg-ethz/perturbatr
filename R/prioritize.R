#' Select hits from an analyzed data set
#'
#' @export
#' @import data.table
#'
#' @param obj  the data.table to be called
#' @param ...  additional parameters
#'  \itemize{
#'  \item{\emph{hit.ratio} }{ the ratio of siRNAs hits such that a gene counts as hit}
#'  \item{\emph{readout.thresh} }{ the mean readout threshold for all the siRNAs}
#'  \item{\emph{p.value.thresh} }{ does a log transformation}
#' }
prioritize <-
function
(
  obj,
  ...
)
{
  UseMethod("prioritize", obj)
}

#' @noRd
#' @export
#' @import data.table
prioritize.svd.analysed.tt <-
function
(
  obj,
  ...
)
{
  res <- .select.hits.tt(obj, ...)
  class(res) <- c("svd.prioritized.tt", "svd.prioritized", class(res))
  invisible(res)
}

#' @noRd
#' @export
#' @import data.table
prioritize.svd.analysed.hyper <-
function
(
  obj,
  ...
)
{
  res <- .select.hits.hyper(obj, ...)
  class(res) <- c("svd.prioritized.hyper", "svd.prioritized", class(res))
  invisible(res)
}

#' @noRd
#' @export
#' @import data.table
prioritize.svd.analysed.pmm <-
function
(
  obj,
  ...
)
{
  res <- .select.hits.pmm(obj, ...)
  class(res) <- c("svd.prioritized.pmm", "svd.prioritized", class(res))
  invisible(res)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr group_by summarize ungroup filter select mutate
.select.hits.tt <-
function
(
  obj,
  ...
)
{
  # TODO here: what do do with multiple sirnas? same hit criterion as in hyper
  params <- list(...)
  hit.rat <- ifelse(hasArg(hit.ratio), params$hit.ratio, 0.5)
  read.thresh <- ifelse(hasArg(readout.thresh), params$readout.thresh, 0.0)
  p.val.thresh <- ifelse(hasArg(p.value.thresh), params$p.value.thresh, 0.05)
  message(paste("Prioritizing on hit.ratio ", hit.rat,
                ", readout threshold " , read.thresh,
                " and p-value threshold ", p.val.thresh, sep=""))
  res <- dplyr::group_by(obj, Virus, Screen, Library,
                         InfectionType, ReadoutType,
                         Design, Cell,
                         GeneSymbol, Entrez,
                         Plate, RowIdx, ColIdx, siRNAIDs) %>%
    dplyr::mutate(Hit=(Pval <= p.val.thresh & abs(Readout) >= read.thresh)) %>%
    ungroup %>%
    dplyr::group_by(Virus, Screen, Library,
                    ReadoutType, InfectionType,
                    Design, Cell,
                    GeneSymbol, Entrez) %>%
    dplyr::summarize(HitRatio = (sum(Hit)/n()),
                     MeanEffect=mean(Readout),
                     MeanPvalue=mean(Pval)) %>%
    ungroup %>%
    dplyr::filter(HitRatio >= hit.rat)
  res
}

#' @noRd
#' @import data.table
#' @importFrom dplyr group_by summarize ungroup filter select
.select.hits.hyper <-
function
(
  obj,
  ...
)
{
  params <- list(...)
  hit.rat      <- ifelse(hasArg(hit.ratio), params$hit.ratio, 0.5)
  message(paste("Prioritizing on hit.ratio ", hit.rat, sep=""))
  res <- dplyr::group_by(obj, Virus, Screen, Library,
                         InfectionType, ReadoutType,
                         Design, Cell,
                         GeneSymbol, Entrez) %>%
    dplyr::summarize(HitRatio        = (base::sum(Hit == TRUE, na.rm=T)/n()),
                     MeanEffect          = base::mean(Readout, na.rm=T),
                     MeanPvalue          = base::mean(Pval)) %>%
    ungroup %>%
    dplyr::filter(HitRatio >= hit.rat)
  res
}

#' @noRd
#' @import data.table
#' @importFrom dplyr filter select group_by mutate summarize full_join
.select.hits.pmm <-
function
(
  obj,
  ...
)
{
  params <- list(...)
  fdr.threshold <- ifelse(hasArg(fdr.threshold), params$fdr.threshold, 0.2)
  message(paste("Prioritizing on fdr.threshold ", fdr.threshold, sep=""))
  gpes <- obj$gene.pathogen.matrix
  gene.pathogen.results <- gpes %>%
    dplyr::filter(FDR <= th)
  gene.results   <- gpes %>%
    dplyr::select(GeneSymbol, FDR) %>%
    dplyr::group_by(GeneSymbol) %>%
    dplyr::summarize(FDR=base::min(FDR)) %>%
    ungroup %>%
    dplyr::filter(FDR <= th)
  gene.effect.results <- dplyr::full_join(gene.results,
                                          obj$gene.effects,
                                          by="GeneSymbol")
  res <- list(gene.pathogen.results=gene.pathogen.results,
              gene.effect.results=gene.effect.results)
  res
}
