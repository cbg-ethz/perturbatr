#' @noRd
#' @import data.table
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
.loess <- function(obj)
{
  message("Calculating LOESS model!")
  loessed.data <-
    dplyr::group_by(obj, Virus, Screen, Library,
                    ReadoutType, InfectionType, ReadoutClass,
                    Design, Cell,
                    Replicate, Plate) %>%
    dplyr::mutate(Readout = .loess.plate(NumCells, Readout)) %>%
    ungroup
  invisible(loessed.data)
}

#' @noRd
#' @importFrom stats lowess
.loess.plate <- function(n.cells, readout)
{
  good.idxs <- (!is.na(n.cells) & !is.na(readout))
  sorted.n.cells     <- sort.int(n.cells[good.idxs], index.return = T)
  sorted.n.cells.idx <- sort.int(sorted.n.cells$ix,  index.return = T)
  loessed.readout <- rep(NA_real_, length(readout))
  if (all(is.na(readout)) || all(is.na(n.cells)))
  {
    loessed.readout <- readout
  }
  else
  {  tryCatch({
    res <- stats::lowess(n.cells[good.idxs], readout[good.idxs])
    loessed.readout[good.idxs] <-
      readout[good.idxs] - res$y[sorted.n.cells.idx$ix] },
    warning = function(war) { message(paste("WARNING: ", war, "\n")); },
    error = function(err)   { message(paste("ERROR: ", err,  "\n")); }
  )
  }
  if (all(is.na(loessed.readout))) loessed.readout <- readout
  invisible(loessed.readout)
}
