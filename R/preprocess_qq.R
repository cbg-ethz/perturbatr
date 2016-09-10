#' @noRd
#' @import data.table
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
.qq.norm <- function(obj, level=c("replicate", "plate"))
{
  level <- match.arg(level)
  # TODO implement on replicates
  message(paste("Calculating quantile-quantile normalisation on", level))
  # TODO difference between plate/replicates
  ret <-
    dplyr::group_by(obj, Virus, Library, Screen,
                    ReadoutType, InfectionType, ReadoutClass,
                    Design, Cell) %>%
    dplyr::mutate(Readout = .qq.norm.group(Replicate, Readout,
                                           RowIdx, ColIdx)) %>%
    ungroup
  invisible(ret)
}

#' @noRd
#' @import data.table
#' @importFrom limma normalizeBetweenArrays
.qq.norm.group <- function(grps, readouts, row.idxs, col.idxs)
{
  if (length(unique(grps)) == 1)
    stop("Please provide multiple replicates/plates per screen!")
  grps.sizes    <- unname(table(grps))
  if (length(unique(grps.sizes)) != 1)
    stop("Please provide data with equal number of samples in every group.")
  ord <- data.table::data.table(Grp=paste("g", grps, sep=""),
                                Readouts=readouts,
                                Row=row.idxs, Col=col.idxs) %>%
    dplyr::group_by(Grp) %>%
    dplyr::mutate(ord=1:n()) %>%
    ungroup %>%
    tidyr::spread(Grp, Readouts)
  qq.mat <- dplyr::select(ord, -Row, -Col, -ord) %>% as.matrix
  res <- ord
  tryCatch ({
    norm.mat <- limma::normalizeBetweenArrays(qq.mat, method="quantile")
    res <- dplyr::select(ord, Row, Col, ord) %>% cbind(norm.mat) %>%
      as.data.table
  }, warning = function(war) { print(paste("WARNING: ", war)); },
  error = function(err) { print(paste("ERROR: ", err)); }
  )
  ret <- tidyr::gather(res, Grp, Readout, g1:g9) %>%
    as.data.table %>%
    .[, .SD[order(ord)], by="Grp"]
  invisible(ret$Readout)
}
