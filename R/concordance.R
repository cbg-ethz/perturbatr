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
