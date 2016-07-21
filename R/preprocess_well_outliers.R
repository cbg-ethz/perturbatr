#' @noRd
#' @import data.table
.remove.outliers <-
function
(
  obj,
  outlier.wells,
  outlier.well.range
)
{
  if (!is.na(outlier.wells))
    obj <- .rm.outlier.wells(obj, outlier.wells, outlier.well.range)
  obj
}

#' @noRd
#' @import data.table
#' @importFrom dplyr mutate group_by
.rm.outlier.wells <-
function
(
  obj,
  outlier.wells,
  outlier.well.range
)
{
  if (outlier.wells == "quantile")
  {

    probs<- c(.05, .95)
    if (all(!is.na(outlier.well.range)))
      probs <- outlier.well.range
    message(paste("Removing wells based on ", probs[1],
                  "% and ", probs[2], "% quantile of cell number!"))
    obj <- dplyr::group_by(obj, Virus, Screen, Library,
                           InfectionType, ReadoutType, ReadoutClass,
                           Cell, Design) %>%
      dplyr::mutate(Readout=.rm.outlier.quantile(Readout, NumCells, probs)) %>%
      ungroup
  }
  else
    stop("Please give a valid method!!")
  obj
}

#' @noRd
#' @importFrom stats quantile
.rm.outlier.quantile <-
function
(
  read,
  num,
  probs
)
{
  re <- read
  if (!all(is.na(num)))
  {
    quant <- unname(quantile(num, na.rm=T, probs=probs))
    message(paste("\t..removing x<" ,quant[1], " | ", quant[2], "<x wells!", sep=""))
    re[num < quant[1] | quant[2] < num] <- NA_real_
  }
  re
}

