#' @noRd
#' @import data.table
.poc <- function(obj, method, ctrl.gene)
{
  message("Calculating POC..")
  if (!is.na(ctrl.gene)) message(paste("..on gene ", ctrl.gene, "!", sep=""))
  f <- .summarization.method(method)
  ret <-
    dplyr::group_by(obj, Virus, Screen, Library,
                    ReadoutType, InfectionType, ReadoutClass,
                    Design, Cell, Replicate, Plate) %>%
    dplyr::mutate(Readout=.poc.grp(Readout, ctrl.gene, f, Control, GeneSymbol))
  ret
}

#' @noRd
.poc.grp <- function(read, ctrl.gene, f, ctrl, genes)
{
  idxs <- which(ctrl == -1)
  if (!is.na(ctrl.gene)) idxs <- which(ctrl == -1 & ctrl.gene == genes)
  if (length(idxs) == 0) stop("No normalization genes found!")
  ret <- (read / f(read[idxs])) * 100
  ret
}
