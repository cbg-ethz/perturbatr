#' @noRd
#' @importFrom RankAggreg RankAggreg
.rankaggreg <- function(obj, genes)
{
  ranks <- matrix(NA, nrow(obj), ncol(obj))
  ranks <- t(apply(obj, 2, function(el)
  {
      tab <- data.table::data.table(Gene=genes, stat=el) %>%
        .[order(el, decreasing=T)]
      tab$Gene
  }))
  RankAggreg::RankAggreg
}
