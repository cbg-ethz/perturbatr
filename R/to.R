#' Calculcate the PMM-prioritized effect matrices for gene and pathogen-gene matrices
#'
#' @export
#' @import data.table
#' @param obj  the object to calculate the effect matrices for
#' @param ...  additional parameters
effect.matrices <- function(obj, ...)
{
  UseMethod("effect.matrices")
}

#' @noRd
#' @export
#' @import data.table
#' @importFrom dplyr filter select
effect.matrices.svd.prioritized.pmm <- function(obj, ...)
{

  gene.top <- obj$gene.effect.hits %>%
    dplyr::select(GeneSymbol, Effect) %>%
    .[order(-abs(Effect))] %>%
    .[, .SD[1:min(25,.N)]]
  pgs <- obj$fit$gene.pathogen.effects %>%
    dplyr::filter(GeneSymbol %in% gene.top$GeneSymbol) %>%
    dplyr::select(Virus, GeneSymbol, Effect) %>%
    tidyr::spread(Virus, Effect)
  res <- list(gene.effects=gene.top, gene.pathogen.hits=pgs)
  class(res) <- "svd.prioritized.pmm.effect.matrices"
  res
}

#' Get the readout matrix (plus control indexes) as list from and svd.plate object
#'
#' @export
#'
#' @param obj  the object for which you want to have the readout matrix
#' @param ...  additional params
readout.matrix <- function(obj, ...) UseMethod("readout.matrix")

#' @noRd
#' @export
#' @import data.table
#' @importFrom dplyr filter
#' @importFrom dplyr select
readout.matrix.svd.plate <- function(obj, ...)
{
  col.size <- max(obj$ColIdx)
  row.size <- max(obj$RowIdx)
  m <- idx <- genes <- matrix(0, row.size, col.size)
  for (i in 1:row.size)
  {
    row <- dplyr::filter(obj, RowIdx==i)
    for (j in 1:col.size)
    {
      col <- dplyr::filter(row, ColIdx==j)
      if (nrow(col) == 0) {
        m[i, j] <- NA_real_
        idx[i, j] <- 0
        genes[i,j] <- NA_character_
      }  else {
        m[i, j]   <- dplyr::select(col, Readout) %>% unlist
        idx[i, j] <- dplyr::select(col, Control)  %>% unlist
        genes[i, j] <- dplyr::select(col, GeneSymbol) %>% unlist
      }
    }
  }
  res <- list(readout=m, idx=idx, genes=genes)
  res
}

#' Get the replicates of a data-set
#'
#' @export
#'
#' @param obj  an object for which the replicates should be retrieved
#' @param ...  additional params
replicates <- function(obj, ...) UseMethod("replicates")

#' @noRd
#' @export
#'
#' @import data.table
#' @importFrom dplyr filter
replicates.svd.raw <- function(obj, ...)
{
  obj <- dplyr::filter(obj, ReadoutClass=="Readout") %>% as.data.table
  replicates.svd.data(obj, ...)
}

#' @noRd
#' @export
#'
#' @import data.table
#' @importFrom dplyr filter
replicates.svd.data <- function(obj, ...)
{
  g <-
    dplyr::group_indices(obj, Virus, Screen, Replicate,
                         Design, Library,
                         ReadoutType, ScreenType) %>%
    as.data.table
  rep.frame <- obj
  rep.frame$grp <- g
  grps <- unique(rep.frame$grp)
  ret <- lapply(grps, function(i)
  {
    pl <- data.table::as.data.table(dplyr::filter(rep.frame, grp==i))
    class(pl) <- c("svd.replicate", class(pl))
    pl
  })
  class(ret) <- "svd.replicates"
  invisible(ret)
}

#' Get the plates of a data-set
#'
#' @export
#'
#' @param obj  an object for which the plates are going to be retrieved
#' @param ...  additional params
plates <- function(obj, ...)
{
  UseMethod("plates", obj)
}

#' @noRd
#' @export
#'
#' @import data.table
#' @importFrom dplyr filter
plates.svd.raw <- function(obj, ...)
{
  obj <- dplyr::filter(obj, ReadoutClass=="Readout") %>%
    as.data.table
  plates.svd.data(obj, ...)
}

#' @noRd
#' @export
#'
#' @import data.table
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
#' @importFrom dplyr group_indices
plates.svd.data <- function(obj, ...)
{
  g <-
    dplyr::group_indices(obj, Virus, Screen, Replicate, Plate,
                         ReadoutType, ScreenType) %>%
    as.data.table
  plate.frame <- obj
  plate.frame$grp <- g
  grps <- unique(plate.frame$grp)
  plates <- lapply(grps, function(i)
  {
    pl <- data.table::as.data.table(dplyr::filter(plate.frame, grp==i))
    class(pl) <- c("svd.plate", class(pl))
    pl
  })
  class(plates) <- "svd.plates"
  invisible(plates)
}


#' Calculate the correlation between two data-sets
#'
#' @export
#' @import data.table
#'
#' @param x  a random object
#' @param y  a random object of the same type (optional)
#' @param method  a character string indicating which correlation coefficient
#' (\emph{pearson}, \emph{kendall} or \emph{spearman})
#' @param ...  additional parameters
correlation <- function(x, y, method = c("pearson", "kendall", "spearman"), ...)
{
  UseMethod("correlation")
}

#' @export
#' @import data.table
#' @importFrom dplyr filter
#' @importFrom stats cor
#' @method correlation svd.replicates
correlation.svd.replicates <- function(x, y,
                                       method = c("pearson", "kendall", "spearman"),
                                       ...)
{
  method <- match.arg(method)
  len <- length(x)
  gri <- expand.grid(1:len, 1:len)
  ret <- do.call(
    "rbind",
    apply(gri, 1, function (gr)
    {
      val1 <- unname(gr[1])
      val2 <- unname(gr[2])
      x1 <- unname(x[[val1]]$Readout)
      x2 <- unname(x[[val2]]$Readout)
      idx <- which(!is.na(x1) & !is.na(x2))
      df <- data.table::data.table(
        Var1=val1, Var2=val2,
        Cor=stats::cor(x1[idx], x2[idx], method=method))
      df
    })
  )
  ret <- data.table::dcast(ret, Var1 ~ Var2, value.var="Cor") %>%
    dplyr::select(-Var1) %>%
    as.matrix(ret)
  colnames(ret) <- rownames(ret) <- paste("Replicate", 1:len, sep="_")
  ret
}
