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
  obj <- dplyr::filter(obj, ReadoutClass=="Readout") %>%
    as.data.table
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
                         ReadoutType, InfectionType) %>%
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
plates <- function(obj, ...) UseMethod("plates", obj)

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
                         ReadoutType, InfectionType) %>%
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

#' @noRd
remove.suffix <-
function(s)
{
  if (length(grep("\\.[[:alpha:]]+$", s, perl=T)) != 0)
  {
    s <- sub("\\.[[:alpha:]]+$", "", s)
  }
  s
}

#' @noRd
is.whole <- function(a) is.numeric(a) && floor(a) == a

#' @noRd
concat <-
function(arr, sep, col)
{
  if (length(arr) == 1) return(arr)
  else paste(arr, sep=sep, col=col)
}

#' @noRd
#' @importFrom stats median
.summarization.method <- function(summ.method)
{
    f <- switch(as.character(summ.method),
                "mean"=base::mean,
                "median"=stats::median,
                "min"=base::min,
                `NA`=NA,
                stop("wrong method given"))
    f
}

#' Calculate the concordance between the vectorial elements of a list.
#'
#' @export
#'
#' @docType methods
#' @rdname concordance-methods
#'
#' @param obj  a list of vectors of the same type
#' @param ...  additional parameters
setGeneric(
  "concordance",
  function(obj, ...) standardGeneric("concordance")
)

#' @rdname concordance-methods
#' @aliases concordance,list-method
setMethod(
  "concordance",
  c(obj="list"),
  function(obj, ...)
  {
    if (length(obj) < 0) stop("Please provide some arguments!")
    if (is.null(names(obj)) |
        any(is.null(names(obj))) |
        any(names(obj) == "")) stop("Please provide names for list items!")
    classes <- unname(sapply(obj, class))
    if (any(classes == "list"))
      stop("Please do not provide lists as elements!")
    if (any(classes != classes[1]))
      stop("Please provide the same class for every list element!")
    invisible(concordance.default(obj, ...))
  }
)

#' @noRd
#'
#' @importFrom data.table melt
concordance.default <- function(obj, ...)
{
  oo.m <- matrix(0, length(obj), length(obj))
  jac.m <- matrix(0, length(obj), length(obj))
  for(i in 1:length(obj))
  {
    el1 <- obj[[i]]
    for (j in 1:length(obj))
    {
      el2 <- obj[[j]]
      jac.m[i, j] <- length(intersect(el1, el2)) /
        length(union(el1, el2))
      oo.m[i, j]  <- length(intersect(el1, el2))  /
        min(length(el1), length(el2))
    }
  }
  oo.df  <- data.table::melt(oo.m)
  jac.df <- data.table::melt(jac.m)
  coord  <-
    list(jaccard=jac.df,
         overlap=oo.df,
         jaccard.matrix=jac.m,
         overlap.matrix=oo.m,
         obj=obj,
         names=names(obj))
  class(coord) <- "svd.concordance"
  invisible(coord)
}

#' Drops columns that are not needed
#'
#' @export
#'
#' @param obj  object that should be dropped
#' @param ... additional params
drop <- function(obj, ...) UseMethod("drop")

#' @noRd
#' @export
drop.svd.data <- function(obj, ...)
{
  cls <- colnames(obj)
  ret <- obj
  if ("Viability" %in% cls) ret <- dplyr::select(ret, -Viability)
  ret <- dplyr::filter(ret, !is.na(GeneSymbol), GeneSymbol != "buffer")
  class(ret) <- class(obj)
  ret
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
cor <- function(x, y, method = c("pearson", "kendall", "spearman"), ...)
  UseMethod("cor")

#' @noRd
#' @export
#' @import data.table
#' @importFrom dplyr filter
#' @importFrom stats cor
cor.svd.replicates <- function(x, y,
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


#' Calculcate the I-dont-know-how-to-call-them
#'
#' @export
#' @import data.table
#' @param obj  the object to calculate the ffect matrices for
#' @param ...  additional parameters
effect.matrices <- function(obj, ...) UseMethod("effect.matrices")

#' @noRd
#' @export
#' @import data.table
#' @importFrom dplyr filter select
effect.matrices.svd.prioritized.pmm <- function(obj, ...)
{
  gene.hits <- data.table::as.data.table(obj$res$all.virus.results) %>%
    .[order(abs(Effect), decreasing=T)]
  pathogen.gene.hits <-
    base::do.call("rbind",
                  base::lapply(obj$res$single.virus.results , function(i) i))
  pathogen.gene.hits$Virus <-
    base::unname(unlist(sapply(rownames(pathogen.gene.hits),
                               function(e) sub(".[[:digit:]]+", "", e))))
  cpgs <- data.table::as.data.table(obj$model$gene.pathogen.effects)

  cpg.mat <- dplyr::filter(cpgs, GeneID %in% gene.hits$GeneSymbol) %>%
    dplyr::select(GeneID, grep("GenePathogenEffect", colnames(cpgs)))
  base::colnames(cpg.mat) <- c("GeneSymbol", "CHIKV", "DENV", "HCV", "SARS")

  fdr.mat <- dplyr::filter(cpgs, GeneID %in% gene.hits$GeneSymbol) %>%
    dplyr::select(GeneID, grep("fdr", colnames(cpgs)))
  base::colnames(fdr.mat) <- c("GeneSymbol", "CHIKV", "DENV", "HCV", "SARS")

  res        <- base::list(cpg.mat=cpg.mat, fdr.mat=fdr.mat)
  class(res) <- "svd.pmm.single.gene.matrices"
  res
}

#' @noRd
#' @import igraph
#' @importFrom igraph graph.data.frame
#' @importFrom utils read.csv
.read.graph <- function(pth)
{
  tab <- utils::read.csv(pth, sep="\t", header=T)
  gra <- igraph::graph.data.frame(read.csv(pth, sep="\t", header=T),
                                  directed=F)
  if (ncol(tab) == 3)
  {
    if (is.null(E(gra)$weight))
      stop("Third column sould be 'weight'")
  }
  gra
}

#' @noRd
#' @import Matrix
#' @importFrom assertthat assert_that
.stoch.col.norm <- function(m)
{
    col.sums  <- Matrix::colSums(m)
    zero.cols <- Matrix::which(col.sums == 0)
    if (length(zero.cols) != 0)
    {
      m[,zero.cols]       <- 1
      col.sums[zero.cols] <- nrow(m)
    }
    ret <- sweep(m, 2L, col.sums, "/", check.margin = FALSE)
    ret.col.sums <- unname(Matrix::colSums(ret))
    assertthat::assert_that(all(ret.col.sums <= 1.001))
    assertthat::assert_that(all(ret.col.sums >= 0.999))
    assertthat::assert_that(all(ret >= 0) & all(!is.nan(ret@x)))
    invisible(ret)
}
