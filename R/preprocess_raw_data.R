#' @noRd
#' @import data.table
#' @importFrom tidyr spread
.normalize <-
function
(
  obj,
  normalize,
  normalize.viability,
  ...
)
{
  res <-
    .normalize.within(obj, normalize, normalize.viability, ...) %>%
    tidyr::spread(ReadoutClass, Readout)
  res
}

#' @noRd
#' @import data.table
#' @importFrom dplyr filter
.normalize.within <-
function
(
  obj,
  normalize,
  normalize.viability,
  ...
)
{
  if (!("Readout" %in% colnames(obj))) stop("No 'Readout' column given!")
  read.dat <- dplyr::filter(obj, ReadoutClass == "Readout")
  via.dat <- dplyr::filter(obj, ReadoutClass == "Viability")
  # pars arguments
  params <- list(...)
  z.score.level <- ifelse(hasArg(z.score.level),
                          params$z.score.level, "plate")
  z.score.ctrl  <- ifelse(hasArg(z.score.ctrl),
                          params$z.score.ctrl, NA_character_)
  poc.ctrl      <- ifelse(hasArg(poc.ctrl),
                          params$poc.ctrl, NA_character_)
  method        <- ifelse(hasArg(method),
                          params$method, "mean")
  background.column <- ifelse(hasArg(background.column),
                              params$background.column,
                              NA)
  background.row <- ifelse(hasArg(background.row), params$background.row, NA)
  # normalize the readout
  message("Normalizing readout")
  read.dat <- .do.normalize.within(obj=read.dat,
                                  normalize=normalize,
                                  method=method,
                                  poc.ctrl=poc.ctrl,
                                  z.score.level=z.score.level,
                                  z.score.ctrl=z.score.ctrl,
                                  background.column=background.column,
                                  background.row=background.row)
  # normalize the viability
  if (normalize.viability)
  {
    message("Normalizing viability")
    via.dat <- .do.normalize.within(obj=via.dat,
                                     normalize=normalize,
                                     method=method,
                                     poc.ctrl=poc.ctrl,
                                     z.score.level=z.score.level,
                                     z.score.ctrl=z.score.ctrl,
                                     background.column=background.column,
                                     background.row=background.row)
  }
  invisible(data.table::rbindlist(list(read.dat, via.dat)))
}

#' @noRd
#' @import data.table
.do.normalize.within <-
function
(
  obj,
  normalize,
  method,
  poc.ctrl,
  z.score.level,
  z.score.ctrl,
  background.column,
  background.row
)
{
  for (norm in normalize)
  {
    obj <- switch(norm,
                  "log"=.log.norm(obj=obj),
                  "poc"=.poc(obj=obj, method=method, poc.ctrl=poc.ctrl),
                  "z.score"=.z.score(obj=obj, method="default",
                                      level=z.score.level,
                                      ctrl=z.score.ctrl),
                  "robust-z.score"=.z.score(obj=obj, method="robust",
                                            level=z.score.level,
                                            ctrl=z.score.ctrl),
                  "loess"=.loess(obj=obj),
                  "b.score"=.b.score(obj=obj),
                  "background"=.background.correct(
                    obj=obj, background.column=background.column,
                    background.row=background.row, method=method),
                  "qq"=.qq.norm(obj=obj),
                  stop("Give a default normalization method!"))
  }
  obj
}

#' @noRd
#' @import data.table
#' @importFrom dplyr mutate group_by select
.background.correct <-
function
(
  obj,
  background.column,
  background.row,
  method
)
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
.substract.background <-
function
(
  readout,
  f,
  background
)
{
  ret  <- readout - f(readout[background])
  ret
}


#' @noRd
#' @import data.table
.poc <-
function
(
  obj,
  method,
  ctrl.gene
)
{
  message("Calculating POC..")
  if (!is.na(ctrl.gene)) message(paste("..on gene ", ctrl.gene, "!", sep=""))
  f <- .summarization.method(method)
  ret <-
    dplyr::group_by(obj, Virus, Screen, Library,
                  ReadoutType, InfectionType, ReadoutClass,
                  Design, Cell,
                  Replicate, Plate) %>%
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

#' @noRd
#' @import data.table
#' @importFrom dplyr mutate
.log.norm <-
function
(
 obj
)
{
  message("Calculating log!")
  obj <- dplyr::mutate(obj, Readout = log(Readout  + 0.00001))
  invisible(obj)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
.z.score <-
function
(
 obj,
 method=c("default", "robust"),
 level=c("plate", "control"),
 ctrl
)
{
  method <- match.arg(method)
  level  <- match.arg(level)
  message(paste("Calculating z-scores on", level))
  if (!is.na(ctrl)) message(paste("...normalizing on", ctrl))
  z.score.data <-
    dplyr::group_by(obj, Virus, Screen, Library,
                    ReadoutType, InfectionType, ReadoutClass,
                    Design, Cell,
                    Replicate, Plate)
  if (level == "plate")
  {
    z.score.data <-
      dplyr::mutate(z.score.data, Readout = .z.score.plate(Readout, method))
  }
  else if (level == "control")
  {
    check <- dplyr::filter(z.score.data, Control == -1) %>%
      dplyr::mutate(n=n())
    if (length(unique(check$n))  > 1)
      warning(paste("Found unequal number of negative controls on plates:",
                    paste(unique(check$n), collapse=", ")))
    z.score.data <-
      dplyr::mutate(z.score.data,
                    Readout = .z.score.control(Readout, method,
                                               Control, GeneSymbol, ctrl))
  }
  dplyr::ungroup(z.score.data)
  invisible(z.score.data)
}

#' @noRd
#'
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom stats mad
.z.score.plate <-
function
(
  obj,
  method
)
{
  val <- switch(method,
                "default"=((obj - base::mean(obj, na.rm=T)) /
                             stats::sd(obj, na.rm=T)),
                "robust" =((obj - stats::median(obj, na.rm=T)) /
                             stats::mad(obj, na.rm=T)),
                stop("No correct method given!"))
  invisible(val)
}

#' @noRd
#'
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom stats mad
.z.score.control <-
function
(
  re,
  method,
  control,
  genes,
  ctrl.gene
)
{
  genes <- tolower(genes)
  ctrl.gene <- tolower(ctrl.gene)
  cont.idx <- which(control == -1)
  if (!is.na(ctrl.gene)) {
      cont.idx <- which(control == -1 & genes == ctrl.gene)
      if (length(cont.idx) == 0) stop("No controls found for criteria!")
  }
  if (is.na(cont.idx[1]))
  {
    warning("There are no negative controls on your plate!
            Stopping control-wise normalization and using plate-wise instead.")
    val <- .z.score.plate(re, method)
  }
  else
  {
    if (length(cont.idx) < 3)
      warning(paste("You are normalizing with z-scores and only using",
                    length(cont.idx),"values!"))
    val <- switch(method,
                  "default"=((re - base::mean(re[cont.idx], na.rm=T)) /
                               stats::sd(re[cont.idx], na.rm=T)),
                  "robust" =((re - stats::median(re[cont.idx], na.rm=T)) /
                               stats::mad(re[cont.idx], na.rm=T)),
                  stop("No correct method given!"))
  }
  val
}

#' @noRd
#' @import data.table
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
.loess <-
function
(
 obj
)
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
.loess.plate <-
function
(
 n.cells,
 readout
)
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

#' @noRd
#' @import data.table
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
.b.score <-
function
(
 obj
)
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
#'
#' @importFrom stats medpolish
#' @importFrom stats median
#' @importFrom stats mad
.bscore.plate <-
function
(
  rows,
  cols,
  readouts
)
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

#' @noRd
#' @import data.table
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
.qq.norm <-
function
(
 obj,
 level=c("replicate", "plate")
)
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
.qq.norm.group <-
function
(
 grps,
 readouts,
 row.idxs,
 col.idxs
)
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

