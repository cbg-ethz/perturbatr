#' @noRd
#' @importFrom RankAggreg RankAggreg
#' @param obj  a n x k matrix where n is the number of lists and k the number of elements in each list
#' @param labels  the names of the k elements
#' @param k  how many elements should be considered
.rankaggreg <- function(obj, labels, k)
{
  ranks <- matrix(NA, nrow(obj), ncol(obj))
  ranks <- t(apply(obj, 2, function(el)
  {
      tab <- data.table::data.table(lab=labels, stat=el) %>%
        .[order(el, decreasing=T)]
      tab$lab
  }))
  ret <- RankAggreg::RankAggreg(ranks[,1:100],
                                method="CE",distance="Kendall", k=k)
  ret
}
