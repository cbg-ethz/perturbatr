#' Fit an LMM to the data and calculate local false discovery rates.
#'
#' @export
#'
#' @import data.table
#'
#' @param obj  an svd.data object
#' @param th  the threshold for the local false discovery rate
#' @param drop  boolean flag if all entries should be dropped that are not found in every virus
#' @param ...  additional parameters
#' \itemize{
#'  \item{ignore }{ remove sirnas that have been found less than \code{ignore} times}
#' }
lmm <- function(obj, th=.1, drop=T, ...) UseMethod("lmm", obj)

#' @export
#' @noRd
#' @import data.table
lmm.svd.data <-
function
(
  obj,
  th=.1,
  drop=T,
  ...
)
{
  lf  <- .lmm(obj, drop)
  tlf <- .thresh(lf, th)
  res <- list(model=lf, res=tlf)
  class(res) <- c("svd.analysed.pmm","svd.analysed", class(res))
  mat <- .set.gene.matrices(res)
  res$result.matrices <- mat
  invisible(res)
}

#' @noRd
#' @import data.table
#' @import lme4
#' @importFrom dplyr mutate
#' @importFrom dplyr select
.lmm <-
function
(
  obj,
  drop,
  ...
)
{
  params <- base::list(...)
  ignore <- base::ifelse(hasArg(ignore) && is.numeric(params$ignore), params$ignore, 1)
  # init the data table for the LMM
  model.data <- .set.lmm.matrix(obj, drop, ignore)
  # save gene control mappings
  gene.control.map <-
    dplyr::select(model.data, GeneSymbol, Control) %>%
    unique
  gene.control.map$GeneSymbol <- as.character(gene.control.map$GeneSymbol)
  # fit the LMM
  fit.lmm <- lme4::lmer(Readout ~ Virus + (1 | GeneSymbol) + (1 | Virus:GeneSymbol),
                        data = model.data, weights = model.data$Weight,
                        verbose = FALSE)
  random.effects <- lme4::ranef(fit.lmm)
  # create the data table with gene effects
  ag <- data.table::data.table(ag = random.effects[["GeneSymbol"]][,1],
                   GeneSymbol = as.character(rownames(random.effects[["GeneSymbol"]])))
  # create the data.table with pathogen effects
  bcg <- data.table::data.table(bcg = random.effects[["Virus:GeneSymbol"]][,1],
                    GenePathID = as.character(rownames(random.effects[["Virus:GeneSymbol"]]))) %>%
    dplyr::mutate(GeneSymbol = sub("^.+:", "", GenePathID))
  # create the table with gene-pathogen effects
  ccg <- base::merge(bcg, ag, by = "GeneSymbol") %>%
    dplyr::mutate(Virus = sub(":.+$", "", GenePathID), GeneVirusEffect = ag + bcg) %>%
    dplyr::select(-GenePathID, -bcg, -ag)
  # calculate fdrs
  fdrs <- .fdr(ccg)
  # finalize output and return as list
  ret <- base::list(
    gene.effects=data.table::as.data.table(
      base::merge(ag, gene.control.map, by="GeneSymbol")),
    gene.pathogen.effects=data.table::as.data.table(
      base::merge(fdrs$ccg.matrix, gene.control.map, by="GeneSymbol")))
  colnames(ret$gene.effects) <- c("GeneSymbol", "Effect", "Control")

  ret$ag <- ag
  ret$bcg <- bcg
  ret$ccg <- ccg
  ret$model.data <- model.data
  ret$fit <- fit.lmm
  ret$fdrs <- fdrs$fdrs
  invisible(ret)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
.set.lmm.matrix <-
function
(
  obj,
  drop,
  ignore
)
{
  pmm.mat <-
    dplyr::select(obj, Entrez, GeneSymbol, Virus, Readout, Control) %>%
    dplyr::filter(!is.na(GeneSymbol))
  data.table::setDT(pmm.mat)[, Weight := 1]
  data.table::setDT(pmm.mat)[, VG := paste(Virus, GeneSymbol, sep=":")]
  pmm.mat$Entrez     <- as.integer(pmm.mat$Entrez)
  pmm.mat$Virus      <- as.factor(pmm.mat$Virus)
  pmm.mat$VG         <- as.factor(pmm.mat$VG)
  pmm.mat$GeneSymbol <- as.factor(pmm.mat$GeneSymbol)
  pmm.mat$Weight     <- as.integer(pmm.mat$Weight)
  pmm.mat$Control    <- as.integer(pmm.mat$Control)
  pmm.mat <-
    dplyr::filter(pmm.mat, !is.na(Readout)) %>%
    dplyr::group_by(VG) %>%
    dplyr::mutate(cnt=n()) %>%
    ungroup %>%
    dplyr::filter(cnt >= ignore) %>%
    dplyr::select(-cnt)
  if (drop)
  {
    vir.cnt <- length(unique(pmm.mat$Virus))
    pmm.mat <- dplyr::group_by(pmm.mat, GeneSymbol) %>%
      dplyr::mutate(drop=(length(unique(Virus)) != vir.cnt)) %>%
      ungroup %>%
      dplyr::filter(!drop) %>%
      dplyr::select(-drop)
  }
  pmm.mat <- droplevels(pmm.mat)
  invisible(pmm.mat)
}

#' @noRd
#' @importFrom stats reshape
#' @importFrom stats na.omit
.fdr <-
function
(
  obj,
  ...
)
{
  ccg.matrix <-
    stats::reshape(as.data.frame(obj),
                   direction = "wide", timevar = "Virus",
                   idvar = "GeneSymbol")
  fdrs <- list()
  for (i in unique(obj$Virus))
  {
    cl <- paste("GeneVirusEffect", i, sep = ".")
    fr <- paste("fdr", i, sep = ".")
    sub.ccg.matrix <- stats::na.omit(ccg.matrix[, c("GeneSymbol",cl) ])
    locf <- .localfdr(sub.ccg.matrix[, cl])
    fdrs[[i]] <- locf
    sub.ccg.matrix <- cbind(sub.ccg.matrix, fdr = locf$fdr)
    nmiss <- sum(is.na(ccg.matrix[, cl]))
    dmiss <- ccg.matrix[is.na(ccg.matrix[, cl ]), c("GeneSymbol", cl)]
    sub.ccg.matrix <- rbind(sub.ccg.matrix, cbind(dmiss, fdr = rep(NA, nmiss)))
    colnames(sub.ccg.matrix)[colnames(sub.ccg.matrix) == "fdr"] <- fr
    ccg.matrix <- merge(ccg.matrix, sub.ccg.matrix[, c("GeneSymbol", fr)],
                        by = "GeneSymbol")
  }
  list(ccg.matrix=ccg.matrix, fdrs=fdrs)
}

#' This is the implementation of Efron local fdr with some additions regarding the return values
#'
#' @noRd
#' @import locfdr
#' @importFrom splines ns
#' @importFrom stats glm quantile poly lm approx poisson qnorm
.localfdr <-
  function
(
  zz,
  bre = 120,
  df = 7, pct = 0, pct0 = 1/4, nulltype = 1,
  type = 0, mult, mlests, main = " ", sw = 0
)
{
  call = match.call()
  if (length(bre) > 1) {
    lo <- min(bre)
    up <- max(bre)
    bre <- length(bre)
  }
  else {
    if (length(pct) > 1) {
      lo <- pct[1]
      up <- pct[2]
    }
    else {
      if (pct == 0) {
        lo <- min(zz)
        up <- max(zz)
      }
      if (pct < 0) {
        med = median(zz)
        ra = med + (1 - pct) * (range(zz) - med)
        lo = ra[1]
        up = ra[2]
      }
      if (pct > 0) {
        v <- quantile(zz, c(pct, 1 - pct))
        lo <- v[1]
        up <- v[2]
      }
    }
  }
  zzz <- pmax(pmin(zz, up), lo)
  breaks <- seq(lo, up, length = bre)
  zh <- hist(zzz, breaks = breaks, plot = F)
  x <- (breaks[-1] + breaks[-length(breaks)])/2
  yall <- y <- zh$counts
  K <- length(y)
  N <- length(zz)
  if (pct > 0) {
    y[1] <- min(y[1], 1)
    y[K] <- min(y[K], 1)
  }
  if (type == 0) {
    X <- cbind(1, splines::ns(x, df = df))
    f <- stats::glm(y ~ splines::ns(x, df = df), stats::poisson)$fit
  }
  else {
    X <- cbind(1, stats::poly(x, df = df))
    f <- stats::glm(y ~ stats::poly(x, df = df), stats::poisson)$fit
  }
  l <- log(f)
  Fl <- cumsum(f)
  Fr <- cumsum(rev(f))
  D <- (y - f)/(f + 1)^0.5
  D <- sum(D[2:(K - 1)]^2)/(K - 2 - df)
  if (D > 1.5)
    warning(paste("f(z) misfit = ", round(D, 1), ".  Rerun with increased df",
                  sep = ""))
  if (nulltype == 3) {
    fp0 = matrix(NA, 6, 4)
    colnames(fp0) = c("delta", "sigleft", "p0", "sigright")
  }
  else {
    fp0 = matrix(NA, 6, 3)
    colnames(fp0) = c("delta", "sigma", "p0")
  }
  rownames(fp0) = c("thest", "theSD", "mlest", "mleSD", "cmest",  "cmeSD")
  fp0["thest", 1:2] = c(0, 1)
  fp0["theSD", 1:2] = 0
  imax <- seq(l)[l == max(l)][1]
  xmax <- x[imax]
  if (length(pct0) == 1) {
    pctup <- 1 - pct0
    pctlo <- pct0
  }
  else {
    pctlo <- pct0[1]
    pctup <- pct0[2]
  }
  lo0 <- stats::quantile(zz, pctlo)
  hi0 <- stats::quantile(zz, pctup)
  nx <- length(x)
  i0 <- (1:nx)[x > lo0 & x < hi0]
  x0 <- x[i0]
  y0 <- l[i0]
  if (nulltype == 3) {
    X00 <- cbind((x0 - xmax)^2, pmax(x0 - xmax, 0)^2)
  }
  else {
    X00 <- cbind(x0 - xmax, (x0 - xmax)^2)
  }
  lr <- stats::lm(y0 ~ X00)
  co <- lr$coef
  if (nulltype == 3) {
    cmerror = I(is.na(co[3]) | is.na(co[2]))
    if (!cmerror)
      cmerror = I(co[2] >= 0 | co[2] + co[3] >= 0)
  }
  else {
    cmerror = is.na(co[3])
    if (!cmerror)
      cmerror = I(co[3] >= 0)
  }
  if (cmerror) {
    if (nulltype == 3)
      stop("CM estimation failed.  Rerun with nulltype = 1 or 2.")
    else if (nulltype == 2)
      stop("CM estimation failed.  Rerun with nulltype = 1.")
    else {
      X0 <- cbind(1, x - xmax, (x - xmax)^2)
      warning("CM estimation failed, middle of histogram non-normal")
    }
  }
  else {
    if (nulltype == 3) {
      X0 <- cbind(1, (x - xmax)^2, pmax(x - xmax, 0)^2)
      sigs <- 1/sqrt(-2 * (c(co[2], co[2] + co[3])))
      fp0["cmest", c(1, 2, 4)] <- c(xmax, sigs)
    }
    else {
      X0 <- cbind(1, x - xmax, (x - xmax)^2)
      xmaxx <- -co[2]/(2 * co[3]) + xmax
      sighat <- 1/sqrt(-2 * co[3])
      fp0["cmest", 1:2] <- c(xmaxx, sighat)
    }
    l0 <- as.vector(X0 %*% co)
    f0 <- exp(l0)
    p0 <- sum(f0)/sum(f)
    f0 <- f0/p0
    fp0["cmest", 3] <- p0
  }
  b = 4.3 * exp(-0.26 * log(N, 10))
  if (missing(mlests)) {
    med = median(zz)
    sc = diff(quantile(zz)[c(2, 4)])/(2 * stats::qnorm(0.75))
    mlests = locfdr:::locmle(zz, xlim = c(med, b * sc))
    if (N > 5e+05) {
      warning("length(zz) > 500,000: For ML estimation, a wider interval than optimal was used.  To use the optimal interval, rerun with mlests = c(",
              mlests[1], ", ", b * mlests[2], ").\n", sep = "")
      mlests = locfdr:::locmle(zz, xlim = c(med, sc))
    }
  }
  if (!is.na(mlests[1])) {
    if (N > 5e+05)
      b = 1
    if (nulltype == 1) {
      Cov.in = list(x = x, X = X, f = f, sw = sw)
      ml.out = locfdr:::locmle(zz, xlim = c(mlests[1], b * mlests[2]),
                               d = mlests[1], s = mlests[2], Cov.in = Cov.in)
      mlests = ml.out$mle
    }
    else mlests = locfdr:::locmle(zz, xlim = c(mlests[1], b * mlests[2]),
                         d = mlests[1], s = mlests[2])
    fp0["mlest", 1:3] = mlests[1:3]
    fp0["mleSD", 1:3] = mlests[4:6]
  }
  if (sum(is.na(fp0[c(3, 5), 1:2])) == 0 & nulltype > 1)
    if (abs(fp0["cmest", 1] - mlests[1]) > 0.05 | abs(log(fp0["cmest",
                                                              2]/mlests[2])) > 0.05)
      warning("Discrepancy between central matching and maximum likelihood estimates.\nConsider rerunning with nulltype = 1")
  if (is.na(mlests[1])) {
    if (nulltype == 1) {
      if (is.na(fp0["cmest", 1]))
        stop("CM and ML Estimation failed, middle of histogram non-normal")
      else stop("ML estimation failed.  Rerun with nulltype=2")
    }
    else warning("ML Estimation failed")
  }
  if (nulltype < 2) {
    delhat = xmax = xmaxx = mlests[1]
    sighat = mlests[2]
    p0 = mlests[3]
    f0 = stats::dnorm(x, delhat, sighat)
    f0 = (sum(f) * f0)/sum(f0)
  }
  fdr = pmin((p0 * f0)/f, 1)
  f00 <- exp(-x^2/2)
  f00 <- (f00 * sum(f))/sum(f00)
  p0theo <- sum(f[i0])/sum(f00[i0])
  fp0["thest", 3] = p0theo
  fdr0 <- pmin((p0theo * f00)/f, 1)
  f0p <- p0 * f0
  if (nulltype == 0)
    f0p <- p0theo * f00
  F0l <- cumsum(f0p)
  F0r <- cumsum(rev(f0p))
  Fdrl <- F0l/Fl
  Fdrr <- rev(F0r/Fr)
  Int <- (1 - fdr) * f * (fdr < 0.9)
  if (sum(x <= xmax & fdr == 1) > 0)
    xxlo <- min(x[x <= xmax & fdr == 1])
  else xxlo = xmax
  if (sum(x >= xmax & fdr == 1) > 0)
    xxhi <- max(x[x >= xmax & fdr == 1])
  else xxhi = xmax
  if (sum(x >= xxlo & x <= xxhi) > 0)
    fdr[x >= xxlo & x <= xxhi] <- 1
  if (sum(x <= xmax & fdr0 == 1) > 0)
    xxlo <- min(x[x <= xmax & fdr0 == 1])
  else xxlo = xmax
  if (sum(x >= xmax & fdr0 == 1) > 0)
    xxhi <- max(x[x >= xmax & fdr0 == 1])
  else xxhi = xmax
  if (sum(x >= xxlo & x <= xxhi) > 0)
    fdr0[x >= xxlo & x <= xxhi] <- 1
  if (nulltype == 1) {
    fdr[x >= mlests[1] - mlests[2] & x <= mlests[1] + mlests[2]] = 1
    fdr0[x >= mlests[1] - mlests[2] & x <= mlests[1] + mlests[2]] = 1
  }
  p1 <- sum((1 - fdr) * f)/N
  p1theo <- sum((1 - fdr0) * f)/N
  fall <- f + (yall - y)
  Efdr <- sum((1 - fdr) * fdr * fall)/sum((1 - fdr) * fall)
  Efdrtheo <- sum((1 - fdr0) * fdr0 * fall)/sum((1 - fdr0) *
                                                  fall)
  iup <- (1:K)[x >= xmax]
  ido <- (1:K)[x <= xmax]
  Eleft <- sum((1 - fdr[ido]) * fdr[ido] * fall[ido])/sum((1 -
                                                             fdr[ido]) * fall[ido])
  Eleft0 <- sum((1 - fdr0[ido]) * fdr0[ido] * fall[ido])/sum((1 -
                                                                fdr0[ido]) * fall[ido])
  Eright <- sum((1 - fdr[iup]) * fdr[iup] * fall[iup])/sum((1 -
                                                              fdr[iup]) * fall[iup])
  Eright0 <- sum((1 - fdr0[iup]) * fdr0[iup] * fall[iup])/sum((1 -
                                                                 fdr0[iup]) * fall[iup])
  Efdr <- c(Efdr, Eleft, Eright, Efdrtheo, Eleft0, Eright0)
  Efdr[which(is.na(Efdr))] = 1
  names(Efdr) <- c("Efdr", "Eleft", "Eright", "Efdrtheo", "Eleft0",
                   "Eright0")
  if (nulltype == 0)
    f1 <- (1 - fdr0) * fall
  else f1 <- (1 - fdr) * fall
  if (!missing(mult)) {
    mul = c(1, mult)
    EE = rep(0, length(mul))
    for (m in 1:length(EE)) {
      xe = sqrt(mul[m]) * x
      f1e = approx(xe, f1, x, rule = 2, ties = mean)$y
      f1e = (f1e * sum(f1))/sum(f1e)
      f0e = f0
      p0e = p0
      if (nulltype == 0) {
        f0e = f00
        p0e = p0theo
      }
      fdre = (p0e * f0e)/(p0e * f0e + f1e)
      EE[m] = sum(f1e * fdre)/sum(f1e)
    }
    EE = EE/EE[1]
    names(EE) = mul
  }
  Cov2.out = locfdr:::loccov2(X, X0, i0, f, fp0["cmest", ], N)
  Cov0.out = locfdr:::loccov2(X, matrix(1, length(x), 1), i0, f, fp0["thest",
                                                                     ], N)
  if (sw == 3) {
    if (nulltype == 0)
      Ilfdr = Cov0.out$Ilfdr
    else if (nulltype == 1)
      Ilfdr = ml.out$Ilfdr
    else if (nulltype == 2)
      Ilfdr = Cov2.out$Ilfdr
    else stop("With sw=3, nulltype must equal 0, 1, or 2.")
    return(Ilfdr)
  }
  if (nulltype == 0)
    Cov = Cov0.out$Cov
  else if (nulltype == 1)
    Cov = ml.out$Cov.lfdr
  else Cov = Cov2.out$Cov
  lfdrse <- diag(Cov)^0.5
  fp0["cmeSD", 1:3] = Cov2.out$stdev[c(2, 3, 1)]
  if (nulltype == 3)
    fp0["cmeSD", 4] = fp0["cmeSD", 2]
  fp0["theSD", 3] = Cov0.out$stdev[1]
  if (sw == 2) {
    if (nulltype == 0) {
      pds = fp0["thest", c(3, 1, 2)]
      stdev = fp0["theSD", c(3, 1, 2)]
      pds. = t(Cov0.out$pds.)
    }
    else if (nulltype == 1) {
      pds = fp0["mlest", c(3, 1, 2)]
      stdev = fp0["mleSD", c(3, 1, 2)]
      pds. = t(ml.out$pds.)
    }
    else if (nulltype == 2) {
      pds = fp0["cmest", c(3, 1, 2)]
      stdev = fp0["cmeSD", c(3, 1, 2)]
      pds. = t(Cov2.out$pds.)
    }
    else stop("With sw=2, nulltype must equal 0, 1, or 2.")
    colnames(pds.) = names(pds) = c("p0", "delhat", "sighat")
    names(stdev) = c("sdp0", "sddelhat", "sdsighat")
    return(list(pds = pds, x = x, f = f, pds. = pds., stdev = stdev))
  }
  p1 <- seq(0.01, 0.99, 0.01)
  cdf1 <- rep(0, 99)
  fd <- fdr
  if (nulltype == 0)
    fd <- fdr0
  for (i in 1:99) cdf1[i] <- sum(f1[fd <= p1[i]])
  cdf1 <- cbind(p1, cdf1/cdf1[99])
  mat <- cbind(x, fdr, Fdrl, Fdrr, f, f0, f00, fdr0, yall,
               lfdrse, f1)
  namat <- c("x", "fdr", "Fdrleft", "Fdrright", "f", "f0",
             "f0theo", "fdrtheo", "counts", "lfdrse", "p1f1")
  if (nulltype == 0)
    namat[c(3, 4, 10)] <- c("Fdrltheo", "Fdrrtheo", "lfdrsetheo")
  dimnames(mat) <- list(NULL, namat)
  z.2 = rep(NA, 2)
  m = order(fd)[nx]
  if (fd[nx] < 0.2) {
    z.2[2] = stats::approx(fd[m:nx], x[m:nx], 0.2, ties = mean)$y
  }
  if (fd[1] < 0.2) {
    z.2[1] = stats::approx(fd[1:m], x[1:m], 0.2, ties = mean)$y
  }
  hist.dat <- list()
  hist.dat$zvalues <- zzz
  hist.dat$yt <- pmax(yall * (1 - fd), 0)
  hist.dat$x <- x
  hist.dat$f <- f
  hist.dat$f0 <- p0 * f0
  if (!is.na(z.2[2])) hist.dat$z.2.2 <- z.2[2]
  if (!is.na(z.2[1])) hist.dat$z.2.1 <- z.2[1]
  if (nulltype == 0) {
    ffdr <- stats::approx(x, fdr0, zz, rule = 2, ties = "ordered")$y
  }
  else ffdr <- stats::approx(x, fdr, zz, rule = 2, ties = "ordered")$y
  vl = list(fdr = ffdr, fp0 = fp0, Efdr = Efdr, cdf1 = cdf1,
            mat = mat, z.2 = z.2)
  if (!missing(mult))
    vl$mult = EE
  vl$call = call
  vl$hist.dat <- hist.dat
  class(vl) <- "svd.pmm.fdr"
  invisible(vl)
}

#' @noRd
#' @import data.table
.thresh <-
function
(
  tlf,
  th,
  ...
)
{
  params       <- list(...)
  m            <- tlf$gene.pathogen.effects
  gene.effects <- tlf$gene.effects
  col.names   <- colnames(m)
  ccg.idx <- grep("ccg", col.names)
  fdr.idx <- grep("fdr", col.names)
  res <- list()
  row.m <- rownames(m[,1])
  for (i in seq_along(ccg.idx))
  {
    vir <- sub("ccg.", "",col.names[ccg.idx[i]])
    fdrs <- m[, fdr.idx[i]]
    idx <- base::which(fdrs <= th & !is.na(fdrs))
    v.r <- m[idx, c(1, ccg.idx[i], fdr.idx[i])]
    colnames(v.r) <- c("Gene", "Effect", "FDR")
    res[[vir]] <- v.r
  }
  narm <- base::ifelse(hasArg(na.rm), params$na.rm, T)
  mat  <- base::as.matrix(m[,fdr.idx])
  rownames(mat) <- row.m
  fdr.mins <-   base::as.matrix(apply(mat, 1, function(e)
  {
    c(min(e, na.rm=narm))
  }), ncol=1)
  res.mat     <- base::data.frame(GeneSymbol=names(fdr.mins[,1]),  FDR=as.vector(fdr.mins[,1]))
  gene.matrix <- base::merge(res.mat, gene.effects,by="GeneSymbol")
  list(all.virus.screen=gene.matrix,
       all.virus.results=gene.matrix[gene.matrix$FDR <= th,],
       single.virus.results=res)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr filter
#' @importFrom dplyr select
.set.gene.matrices <-
  function
(
  obj,
  ...
)
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
