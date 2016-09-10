#' @noRd
#' @import data.table
#' @importFrom dplyr mutate group_by select
.background.correct <- function(obj, background.column, background.row, method)
{
  ret <-
    dplyr::group_by(obj, Virus, Screen, Library,
                    InfectionType, ReadoutType, ReadoutClass,
                    Design, Cell,
                    Replicate, Plate)
  f <- .summarization.method(method)
  if (!is.na(background.row))
  {
    if (!is.numeric(background.row)) stop("Provide numeric index")
    if (background.row > max(ret$RowIdx))
      stop("Please provide a row index that fits the plate!")
    message(paste("Substracting", method,
                  "of row", background.row, " rom every plate!"))
    ret <- dplyr::mutate(ret, bk = (RowIdx == background.row))
  }
  else if (!is.na(background.column))
  {
    if (!is.numeric(background.column))
      stop("Provide numeric index")
    if (background.column > max(ret$ColIdx))
      stop("Please provide a column index that fits the plate!")
    message(paste("Substracting", method,
                  "of column", background.column, "from every plate!"))
    ret <- dplyr::mutate(ret, bk = (ColIdx == background.column))
  }
  else
  {
    message(paste("Substracting", method, "of all NA genes!"))
    ret <- dplyr::mutate(ret, bk = (is.na(GeneSymbol)))
  }
  ret <-
    dplyr::mutate(ret, Readout = .substract.background(Readout, f, bk)) %>%
    ungroup %>%
    dplyr::select(-bk)
  ret
}

#' @noRd
.substract.background <- function(readout, f, background)
{
  ret  <- readout - f(readout[background])
  ret
}
