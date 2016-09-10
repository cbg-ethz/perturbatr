#' @noRd
#' @import data.table
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
.b.score <- function(obj)
{
  message("Calculating b-scores!")
  b.score.data <-
    dplyr::group_by(obj, Virus, Screen, Library,
                    ReadoutType, InfectionType, ReadoutClass,
                    Design, Cell,
                    Replicate, Plate) %>%
    dplyr::mutate(Readout = .bscore.plate(RowIdx, ColIdx, Readout)) %>%
    ungroup
  invisible(b.score.data)
}

#' @noRd
#' @importFrom stats medpolish
#' @importFrom stats median
#' @importFrom stats mad
.bscore.plate <- function(rows, cols, readouts)
{
  nrow <- max(rows)
  ncol <- max(cols)
  fr <- data.frame(rows, cols, readouts, ord=1:length(readouts))
  fr <- fr[order(fr$rows, fr$cols),]
  # this is needed in case the data is not sorted
  readout.mat <- matrix(fr$readouts, nrow, ncol, byrow=T)
  res <- rep(NA_real_, length(readouts))
  tryCatch ({
    med.pol <- stats::medpolish(readout.mat, maxiter=100,
                                na.rm=T, trace.iter=F)$residuals
    res <- as.vector(t(med.pol)) / stats::mad(med.pol, na.rm=T)
  }, warning = function(war) { message(paste("WARNING: ", war)); },
  error = function(err)   { message(paste("ERROR: ", err));   }
  )
  if (all(is.na(res))) res <- fr$readouts
  fr$res <- res
  fr <- fr[order(fr$ord),]
  invisible(fr$res)
}
