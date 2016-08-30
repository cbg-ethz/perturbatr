#' @noRd
#' @export
#' @import data.table
#' @importFrom graphics hist
hist.svd.plate <-
function
(
  x,
  ...
)
{
  .hist(x, ...)
}

#' @noRd
#' @export
#' @import data.table
#' @importFrom dplyr filter
hist.svd.raw <-
function
(
  x,
  ...
)
{
  ret <- dplyr::filter(x, ReadoutClass=="Readout")
  .hist(ret, ...)
}

#' @noRd
#' @export
#' @import data.table
hist.svd.data <-
function
(
  x,
  ...
)
{
 .hist(x, ...)
}

#' @noRd
#' @import ggplot2
#' @import data.table
.hist <-
function
(
  x,
  ...
)
{
  df <- data.frame(x=x$Readout)
  pl <- ggplot2::ggplot(df) +
    ggplot2::geom_histogram(aes(x=x, y=..density..),
                            alpha=.15, position="identity", bins=42) +
    ggplot2::geom_density(aes(x)) +
    theme_bw() +
    theme(text = element_text(size = 14, family = "Helvetica")) +
    xlab("Readout") +
    ylab("Density")
  pl
}


#' @noRd
#' @export
#'
#' @import data.table
plot.svd.raw <-  function(x, y, ...)
{
  x <- dplyr::filter(x, ReadoutClass=="Readout")
  plot.svd.data(x, ...)
}

#' @noRd
#' @export
#' @import ggplot2
#' @import data.table
#'
#' @importFrom dplyr summarize
#' @importFrom dplyr group_by
#' @importFrom tidyr gather
plot.svd.data <-
function
(
  x,
  ...
)
{
  numb.frame <-
    dplyr::group_by(x, Virus, Screen) %>%
    dplyr::summarize(Replicates=length(unique(Replicate)),
                     Genes=length(unique(GeneSymbol))) %>%
    tidyr::gather(Type, Count, Replicates, Genes)

  pl <-
    ggplot2::ggplot(numb.frame, aes(x=Virus, y = Count)) +
    ggplot2::geom_bar(aes(fill=Virus), stat="identity") +
    ggplot2::facet_grid(Type ~ Screen, scales='free_y') +
    ggplot2::scale_fill_brewer(palette="Spectral") +
    ggplot2::geom_text(aes(label = Count, y = Count), size = 3, vjust=-.25) +
    ggplot2::theme_bw()

  pl
}

#' @noRd
#' @export
#' @import ggplot2
#' @import data.table
#' @importFrom RColorBrewer brewer.pal
plot.svd.concordance <-
function
(
  x, y,
  ...
)
{
  params <- list(...)
  type <- ifelse(hasArg(type), params$type, "jaccard")
  size <- ifelse(hasArg(size), params$size, 14)
  mrix <- switch(type,
                 "jaccard"=x$jaccard.matrix,
                 "overlap"=x$overlap.matrix,
                 stop("Wrong type given, give either: jaccard/overlap!"))
  mrix[upper.tri(mrix)] <- NA
  df <- data.table::melt(mrix)
  pl <-
    ggplot2::ggplot(df, aes(factor(Var1), factor(Var2), fill = value)) +
    ggplot2::geom_tile(aes(fill = value), colour = "black") +
    ggplot2::scale_x_discrete(expand = c(0,0), labels=x$names) +
    ggplot2::scale_y_discrete(expand = c(0,0), labels=x$names) +
    ggplot2::scale_fill_gradient2(low = "white",
                                  high = "#3182bd",
                                  na.value="white",
                                  limits = c(0,1),
                                  name="Concordance") +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size = size, family = "Helvetica"),
          aspect.ratio=1,
          axis.title=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank()) +
    ggplot2::ggtitle(ifelse(hasArg("main"), params$main, "")) +
    ggplot2::guides(fill = guide_colorbar(
                                 title.position = "top", title.hjust = 0.5))

  pl
}

#' @noRd
#' @export
#'
#' @import data.table
#' @import ggplot2
#' @importFrom stats cor
plot.svd.replicates <-
  function
(
  x, y,
  ...
)
{
  params <- list(...)
  meth <- ifelse(hasArg(method), params$method, "raw")
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
      if (meth == "qqplot")
      {
        x1 <- sort(x1)
        x2 <- sort(x2)
      }
      idx <- which(!is.na(x1) & !is.na(x2))
      df <- data.table::data.table(
        X=x1, Y=x2,
        RowIdx=val1, ColIdx=val2,
        Cor=stats::cor(x1[idx], x2[idx], method="spearman"))
      df
    })
  )
  ret <- dplyr::filter(ret, RowIdx <= ColIdx)
  ret <- droplevels(ret)
  cors <- dplyr::group_by(ret, RowIdx, ColIdx) %>%
    dplyr::summarize(Cor=mean(Cor)) %>%
    ungroup
  message(paste("Mean correlation:", mean(cors$Cor)))
  cors$Cor <- format(cors$Cor, digits=2, width=1)
  pl <-
    ggplot2::ggplot(ret) +
    ggplot2::geom_point(aes(x=X, y=Y), pch=16, size=1) +
    ggplot2::facet_grid(RowIdx ~ ColIdx) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(ifelse(meth == "qqplot", "QQ-plot", "Scatterplot")) +
    ggplot2::geom_text(data=cors, aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.5,
                      label=paste("r=", Cor, "", sep=""))) +
    ggplot2::theme(axis.title=element_blank())

  pl
}

#' @noRd
#' @export
#'
#' @import data.table
#' @import ggplot2
plot.svd.replicate <-
function
(
  x, y,
  ...
)
{
  params <- list(...)
  if (missing(x) | missing(y)) stop("Please provide two arguments (x and y)!")
  r1 <- x$Readout
  r2 <- y$Readout
  df <- data.frame(X=r1, Y=r2)
  if (meth == "qqplot")
  {
    df$X <- sort(df$X)
    df$Y <- sort(df$Y)
  }
  meth.p <- ifelse(meth == "raw", "Scatterplot", "QQ-plot")
  idx <- which(!is.na(df$X) & !is.na(df$Y))
  cor <- cor(df$X[idx], df$Y[idx], method="spearman")
  pl <-
    ggplot2::ggplot(df) +
    ggplot2::geom_point(aes(x=X, y=Y), pch=16) +
    ggplot2::xlab("Replicate 1") +
    ggplot2::ylab("Replicate 2") +
    ggplot2::ggtitle(paste(meth.p  ," (r=",
                           format(cor, digits=2, width=1), ")", sep=""))  +
    ggplot2::theme_bw()

  pl
}


#' @noRd
#' @export
#'
#' @import data.table
#' @importFrom graphics par
plot.svd.plates <-
function
(
  x, y,
  ...
)
{
  params <- list(...)
  count <- ifelse(hasArg(count), params$count, NA)
  show.controls <- ifelse(hasArg(show.controls) & is.logical(params$show.controls),
                  params$show.controls, T)
  show.gene.names <- ifelse(hasArg(show.gene.names) & is.logical(params$show.gene.names),
                            params$show.gene.names, T)
  graphics::par(ask=T)
  obj.len <- length(x)
  if (is.na(count)) count <- obj.len
  if (obj.len == count) { plates <- seq(obj.len) }
  else { plates <- sort(sample(seq(obj.len), count)) }
  message(paste("Plotting", count, "random plates."), "\n")
  for (i in plates)
  {
    plate <- x[[i]]
    plate.name <- paste("Virus: ", plate$Virus[1],
                        ", screen:", plate$Screen[1],
                        ", replicate:", plate$Replicate[1],
                        ", plate:", plate$Plate[1],
                        ", readout:", plate$ReadoutType[1],
                        ", infection:", plate$InfectionType[1],
                        sep="")
    print(plot.svd.plate(plate, main=plate.name,
                         show.controls=show.controls,
                         show.gene.names=show.gene.names))
  }
  graphics::par(ask=F)
}

#' @noRd
#' @export
#' @import data.table
#' @import ggplot2
#' @importFrom dplyr full_join
plot.svd.plate <-
function
(
  x, y,
  ...
)
{
  params <- list(...)
  show.controls <- ifelse(hasArg(show.controls) & is.logical(params$show.controls),
                          params$show.controls, T)
  show.gene.names <- ifelse(hasArg(show.gene.names) & is.logical(params$show.gene.names),
                          params$show.gene.names, T)
  mat <- readout.matrix(x)
  mat$genes[is.na(mat$genes)] <- ""
  dr  <- data.table::melt(mat$readout)
  di  <- data.table::melt(mat$idx)
  dg  <- data.table::melt(mat$genes)
  df  <- dplyr::full_join(dr, di, by=c("Var1", "Var2"))
  df  <- dplyr::full_join(df, dg, by=c("Var1", "Var2"))
  colnames(df) <- c("Row", "Column", "Readout", "Control", "Gene")
  pl <-
    ggplot2::ggplot(df, aes(x=Column, y=rev(Row)))
  if (show.controls) {
    ctrl <- df$Control
    lwd <- ctrl
    lwd[lwd != 0] <- 1
    pl <- pl + ggplot2::geom_tile(aes(fill=Readout),
                                  color="black",
                                  lwd=lwd)
  } else {
    pl <- pl + ggplot2::geom_tile(aes(fill=Readout), color="black")
  }
  pl <- pl +
    ggplot2::scale_x_discrete(expand = c(0,0),
                              limits=df$Column,
                              name="Column index") +
    ggplot2::scale_y_discrete(expand = c(0,0),
                              limits=df$Row,
                              labels=rev(df$Row),
                              name="Row index") +
    ggplot2::scale_fill_distiller(palette="Spectral", na.value="white") +
    ggplot2::ggtitle(ifelse(hasArg(main), params$main, "")) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size = 12, family = "Helvetica"),
                   aspect.ratio=.75)
  if (show.gene.names) pl <- pl + ggplot2::geom_text(aes(label=df$Gene), size=3)
  pl
}

#' Plot the quality scores
#'
#' @noRd
#' @export
#' @import data.table
#' @import ggplot2
#' @importFrom dplyr select mutate group_indices filter summarize group_by
#' @importFrom tidyr gather
plot.svd.quality <-
function
(
  x, y,
  ...
)
{
  # plot the raw plate values as boxplot
  qual <- x$data
  grps <- dplyr::group_indices(qual, Virus, Screen, Library,
                               InfectionType, ReadoutType,
                               Replicate, Plate)
  df   <- dplyr::mutate(qual, Plate=grps) %>%
    dplyr::select(Virus, Screen, Readout, Plate) %>%
    dplyr::mutate(Plate=as.factor(Plate)) %>%
    dplyr::filter(!is.na(Readout), !is.infinite(Readout))
  pl <-
    ggplot2::ggplot(df, aes(Plate, Readout)) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_grid(Virus ~ Screen) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x=element_blank()) +
    ggplot2::ggtitle("Plate readouts")
  # plot control densities
  df  <- dplyr::select(qual, Virus, Screen, Readout, Control) %>%
    dplyr::filter(!is.na(Readout), !is.infinite(Readout), !is.na(Control), Control != 0) %>%
    dplyr::mutate(Control = as.character(Control))
  data.table::setDT(df)[Control == "-1", Control := "Negative control"]
  data.table::setDT(df)[Control == "1", Control := "Positive control"]
  data.table::setDT(df)[Control == "0", Control := "Normal"]
  pl2 <-
    ggplot2::ggplot(df) +
    ggplot2::geom_histogram(aes(x=Readout, y=..density.., color=Control),
                            alpha=.5, position="identity", bins=42) +
    ggplot2::facet_grid(Virus ~ Screen) +
    ggplot2::theme_bw() +
    ggplot2::ylab("Density") +
    ggplot2::ggtitle("Control densities")
  # plot z factor and ssmd per plate as boxplot
  qual <- x$quality$plate.quality %>% ungroup
  df <- dplyr::select(qual, Virus, Screen, z.fac.control, ssmd) %>%
    tidyr::gather(key, value, z.fac.control, ssmd) %>%
    dplyr::filter(!is.na(value), !is.infinite(value))
  pl3 <-
    ggplot2::ggplot(df) +
    ggplot2::geom_boxplot(aes(key, value), outlier.shape=NA) +
    ggplot2::scale_y_continuous(limits = quantile(df$value, c(0.1, 0.9))) +
    ggplot2::facet_grid(Virus ~ Screen) +
    ggplot2::theme_bw() +
    ggplot2::xlab("Quality metric") +
    ggplot2::ylab("Readout") +
    ggplot2::scale_x_discrete(breaks=c("ssmd", "z.fac.control"),
                     labels=c("SSMD", "Z-factor")) +
    ggplot2::ggtitle("Plate quality measures")
    # plot positive control and negative control values
  qual <- x$data
  df   <- dplyr::mutate(qual, Plate=grps) %>%
    dplyr::select(Readout, Virus, Screen, Plate, Control) %>%
    dplyr::filter(!is.na(Control), Control != 0) %>%
    dplyr::group_by(Virus, Screen, Plate, Control) %>%
    dplyr::summarize(Readout=mean(Readout, na.rm=T)) %>% ungroup %>%
    dplyr::mutate(Plate=as.factor(Plate), Control=as.factor(Control))
  pl4 <-
    ggplot2::ggplot(df) +
    ggplot2::geom_point(aes(x=Plate, y=Readout, color=Control)) +
    ggplot2::theme_bw() +
    ggplot2::geom_line(aes(x=Plate, y=Readout)) +
    ggplot2::facet_grid(Virus ~ Screen) +
    ggplot2::theme(axis.text.x=element_blank()) +
    ggplot2::scale_color_discrete(breaks=c("-1", "1"),
                                  labels=c("Negative control", "Positive control")) +
    ggplot2::ggtitle("Plate controls")

  svd:::.multiplot(pl, pl2, pl3, pl4,  cols=2)
}

#' Plot a barplot showing the number of hits for every virus
#'
#' @noRd
#' @export
#' @import data.table
#' @import ggplot2
#'
#' @importFrom dplyr group_by summarize mutate
plot.svd.prioritized.pmm.single.virus.result <-
function
(
  x, y,
  ...
)
{
  single.res <- do.call("rbind", lapply(x, function(i) i))
  single.res$Virus <- unname(unlist(sapply(rownames(single.res),
                                           function(e) sub(".[[:digit:]]+", "", e))))
  single.res <-
    dplyr::group_by(single.res, Virus, Sign=sign(Effect)) %>%
    dplyr::summarize(cnt=n()) %>%
    ungroup %>%
    dplyr::mutate(Count=cnt*Sign)

  pl <-
    ggplot2::ggplot(single.res, aes(x=Virus, y=Count, fill=Sign)) +
    ggplot2::geom_bar(stat="identity", position="dodge") +
    ggplot2::scale_fill_distiller(palette="Spectral") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="none",
                   text = element_text(size=22),
                   axis.text.y=element_blank()) +
    ggplot2::geom_hline(yintercept=0) +
    ggplot2::ylab("Count hits") +
    ggplot2::geom_text(aes(x=Virus, y=ifelse(Count>0, Count+1, Count-1),
                           label=abs(Count)), size=5, colour="black")

  pl
}

#' Plot a barplot showing the first 25 hits of an integrated PMM screen
#'
#' @noRd
#' @export
#' @import ggplot2
#' @import data.table
#' @importFrom dplyr filter
plot.svd.prioritized.pmm <-
function
(
  x, y,
  ...
)
{
  pl <- .plot.svd.prioritized.pmm(x$gene.effect.hits, main="Gene effects")
  pl2 <- .plot.svd.prioritized.pmm(x$gene.pathogen.effect.hits, main="Gene-virus effects") +
    ggplot2::facet_wrap( ~ Virus, ncol=length(unique(x$gene.pathogen.effect.hits$Virus))/2)
  pl3 <- .multiplot(pl, pl2, cols=2)
  pl3
}

#' @noRd
#' @import data.table
#' @import ggplot2
#' @importFrom dplyr filter
.plot.svd.prioritized.pmm <- function(x, main, ...)
{
  x <- x[order(abs(Effect), decreasing=T), .SD[1:25]] %>%
    dplyr::filter(!is.na(GeneSymbol), !is.na(Effect))
  x.pos.range <- max(abs(x$Effect))
  x.lim  <- c(-x.pos.range, x.pos.range) + c(-x.pos.range, x.pos.range)/5
  pl <-
    ggplot2::ggplot(x) +
    ggplot2::geom_bar(aes(x=GeneSymbol,y=abs(Effect), fill=Effect),
                      stat="identity") +
    ggplot2::scale_fill_distiller(palette="Spectral", limits=x.lim) +
    ylab("Effect") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y=element_blank(),
                   axis.ticks=element_blank()) +
    ggplot2::coord_polar() +
    ggplot2::ggtitle(main)
  pl
}

#' Plot a barplot showing the first 25 hits of the hypergeometric prioritization
#'
#' @noRd
#' @export
#' @import data.table
plot.svd.prioritized.hyper <-
  function
(
  x, y,
  ...
)
{
  .plot.svd.prioritized(x, ...)
}

#' Plot a barplot showing the first 25 hits of the tt prioritization
#'
#' @noRd
#' @export
#' @import data.table
plot.svd.prioritized.tt <-
  function
(
  x, y,
  ...
)
{
  .plot.svd.prioritized(x, ...)
}

#' @noRd
#' @import data.table
#' @import ggplot2
#' @importFrom dplyr group_by summarize mutate filter
.plot.svd.prioritized <-
function
(
  x,
  ...
)
{
  x <- x[order(abs(MeanEffect), decreasing=T), .SD[1:25]] %>%
    dplyr::filter(!is.na(GeneSymbol), !is.na(MeanEffect))
  x.pos.range <- max(abs(x$MeanEffect))
  x.lim  <- c(-x.pos.range, x.pos.range) + c(-x.pos.range, x.pos.range)/5
  pl <-
    ggplot2::ggplot(x) +
    ggplot2::geom_bar(aes(x=GeneSymbol,y=abs(MeanEffect), fill=MeanEffect),
                      stat="identity") +
    ggplot2::scale_fill_distiller(palette="Spectral", limits=x.lim) +
    ylab("Effect") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y=element_blank(),
                   axis.ticks=element_blank()) +
    ggplot2::coord_polar()
  pl
}

#' Plots the the found genes and their cpg values in matrix form
#'
#' @noRd
#' @export
#'
#' @import data.table
#' @import ggplot2
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tidyr gather
plot.svd.prioritized.pmm.single.gene.matrices <-
function
(
  x, y,
  ...
)
{
  params <- list(...)
  size <- ifelse(hasArg(size), params$size, 14)
  dat <- tidyr::gather(x$cpg.mat, Virus, Effect, 2:5)
  LDcolors <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
  or <- rev(x$GeneSymbol)
  pl <-
    ggplot2::ggplot(dat, aes(Virus, GeneSymbol)) +
    ggplot2::geom_tile(aes(fill = Effect), colour=LDcolors[1]) +
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::scale_y_discrete(expand = c(0,0), limits=or, breaks=or, labels=or) +
    ggplot2::scale_fill_gradient2(low=LDcolors[1],
                                  high=LDcolors[11],
                                  na.value=LDcolors[6],
                                  name="Gene-pathogen\neffect") +
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size = 16, family = "Helvetica"),
                   aspect.ratio=2,
                   axis.text.x=element_text(angle=45,  hjust = 1, size=10),
                   axis.text.y=element_text(size=10),
                   axis.title=element_blank(),
                   axis.ticks=element_blank())
  pl
}

#' Plots the densities of the mixture and the null genes and the histogram of the non-null genes
#'
#' @noRd
#' @export
#'
#' @import ggplot2
plot.svd.pmm.fdr <-
function
(
  x, y,
  ...
)
{
  hits <- x$hist.dat
  # the values used for estimation of the densities (f0 and f1)
  zvals <- data.frame(Z=hits$zvalues)
  # frame of estimated non-null hits (yt), density of mixture(f) and density of nulls (f0)
  densities <- data.frame(X=hits$x, y=hits$yt, f=hits$f, f0=hits$f0)
  pl <-
    ggplot2::ggplot(zvals) +
    ggplot2::geom_histogram(aes(Z), bins = 130, alpha = 0.3) +
    ggplot2::geom_bar(data=densities, aes(x=X, y=y), stat = "identity", color="blue") +
    ggplot2::geom_line(data=densities, aes(x=X, y=f0,  colour="p0*f0"), size=1, linetype=1) +
    ggplot2::geom_line(data=densities, aes(x=X, y=f,  colour="f"), size=1, linetype=2) +
    ggplot2::ylab("Frequency") +
    ggplot2::xlab("Gene-pathogen effect") +
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size=14)) +
    ggplot2::scale_colour_manual(name="Density", values=c("green", "blue"))

  pl
}

#' Plots an volcano-plot of the FDR vs readouts
#'
#' @noRd
#' @export
#' @import data.table
#' @importFrom dplyr filter
plot.svd.analysed.pmm <-
function
(
  x, y,
  ...
)
{
  params  <- list(...)
  sig.thresh     <- ifelse(hasArg(sig.thresh), params$fdr.thresh, .2)
  effect.thresh <- ifelse(hasArg(effect.thresh), params$effect.thresh, 0.0)
  log.scale      <- ifelse(hasArg(log.scale), params$log.scale, F)
  gpes <-
    x$gene.pathogen.effects %>%
    dplyr::filter(!is.na(FDR))
  y <- gpes$FDR + 0.0001
  if (log.scale)
  {
    y <- -log(y)
    fdr.thresh <- -log(sig.thresh)
  }
  xl <- ifelse(hasArg(xlab), params$xlab, "Readout")
  yl <- ifelse(log.scale, "-log(fdr)", "Local FDR")
  .plot.svd.analysed(x=gpes$Effect,
                     y=y,
                     ctrls=gpes$Control,
                     genes=paste(gpes$Virus, gpes$GeneSymbol, sep=":"),
                     readout.thresh=effect.thresh,
                     sig.thresh=sig.thresh,
                     xl=xl,
                     yl=yl)
}

#' Plots an volcano-plot of the p-value vs readouts
#'
#' @noRd
#' @export
#' @import data.table
#' @importFrom dplyr filter
plot.svd.analysed <-
function
(
  x, y,
  ...
)
{
  params  <- list(...)
  sig.thresh     <- ifelse(hasArg(sig.thresh), params$sig.thresh, 0.0)
  readout.thresh <- ifelse(hasArg(readout.thresh), params$readout.thresh, 0.0)
  log.scale      <- ifelse(hasArg(log.scale), params$log.scale, F)
  xl <- ifelse(hasArg(xlab), params$xlab, "Readout")
  yl <- ifelse(log.scale, "-log(p-value)", "p-value")
  obj <- x
  y <- obj$Pval + 0.0001
  if (log.scale)
  {
    y <- -log(y)
    sig.thresh <- -log(sig.thresh)
  }
  .plot.svd.analysed(x=obj$Readout,
                     y=y,
                     ctrls=obj$Control,
                     genes=obj$GeneSymbol,
                     readout.thresh=readout.thresh,
                     sig.thresh=sig.thresh,
                     xl=xl,
                     yl=yl)
}

#' @noRd
#' @import ggplot2
.plot.svd.analysed <-
function
(
  x, y,
  ctrls,
  genes,
  readout.thresh,
  sig.thresh,
  xl, yl,
  ...
)
{
  colors <- rep("grey", length(ctrls))
  colors[abs(x) >= readout.thresh & y < sig.thresh] <- "blue"
  colors[which(ctrls == -1)] <- "red"
  colors[which(ctrls == 1)]  <- "green"
  xlim         <- range(x[is.finite(x)])
  ylim         <- range(y[is.finite(y)])
  hits.idx     <- abs(x) >= readout.thresh & y < sig.thresh
  pos.ctrl.idx <- which(ctrls == 1)
  neg.ctrl.idx <- which(ctrls == -1)
  no.hit.idx   <- which(colors == "grey")
  nohits       <- data.frame(Effect=x[no.hit.idx],   sig=y[no.hit.idx],   Gene=genes[no.hit.idx])
  hits         <- data.frame(Effect=x[hits.idx],     sig=y[hits.idx],     Gene=genes[hits.idx])
  pos.ctrls    <- data.frame(Effect=x[pos.ctrl.idx], sig=y[pos.ctrl.idx], Gene=genes[pos.ctrl.idx])
  neg.ctrls    <- data.frame(Effect=x[neg.ctrl.idx], sig=y[neg.ctrl.idx], Gene=genes[neg.ctrl.idx])
  pl <-
    ggplot2::ggplot() +
    ggplot2::xlim(xlim) +
    ggplot2::ylim(ylim) +
    ggplot2::xlab(xl) +
    ggplot2::ylab(yl)
  if (length(nohits$Effect) != 0)
    pl <- pl +
      ggplot2::geom_point(aes(x=nohits$Effect, y=nohits$sig), col="grey",size=.5)
  if (length(hits$Effect) != 0)
    pl <- pl +
      ggplot2::geom_point(aes(x=hits$Effect, y=hits$sig, color="Hit"), size=1.5, alpha=.75)
  if (sig.thresh > 0)
    pl <- pl +
      ggplot2::geom_hline(yintercept=sig.thresh, alpha=.75, linetype="dashed") +
      ggplot2::geom_text(aes(xlim[1], sig.thresh, label = paste("Significance threshold:", sig.thresh)
                             ,vjust = -.25, hjust=0), size = 4)
  if (readout.thresh > 0)
    pl <- pl +
      ggplot2::geom_vline(xintercept=readout.thresh,  alpha=.75, linetype="dotdash") +
      ggplot2::geom_vline(xintercept=-readout.thresh, alpha=.75, linetype="dotdash") +
      ggplot2::geom_text(aes(readout.thresh, ylim[2], label = paste("Effect threshold:", readout.thresh),vjust = 0, hjust=-.01), size = 4)
  if (length(pos.ctrls$Effect) != 0) pl <- pl +
      ggplot2::geom_point(aes(x=pos.ctrls$Effect, y=pos.ctrls$sig, col="Positive control"), size=1.5, alpha=.75)
  if (length(neg.ctrls$Effect) != 0) pl <- pl +
      ggplot2::geom_point(aes(x=neg.ctrls$Effect, y=neg.ctrls$sig, col="Negative control"), size=1.5, alpha=.75)
  if (length(pos.ctrls$Effect) != 0) pl <- pl +
      ggplot2::geom_text(aes(x=pos.ctrls$Effect, y=pos.ctrls$sig, label=pos.ctrls$Gene), hjust=0, vjust=0, check_overlap=TRUE, nudge_x=0.05, size=4)
  if (length(neg.ctrls$Effect) != 0) pl <- pl +
      ggplot2::geom_text(aes(x=neg.ctrls$Effect, y=(neg.ctrls$sig), label=neg.ctrls$Gene), hjust=0, vjust=0, check_overlap=TRUE, nudge_x=0.05, size=4)
  pl <- pl +
    ggplot2::scale_color_manual(name="Points", values=c("blue", "red", "green")) +
    ggplot2::guides(color=guide_legend(title=NULL)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text=element_text(size=12),
                   axis.title=element_text(size=14,face="bold"),
                   legend.text=element_text(size=14))
  pl
}

#' Plots the data used for modelling of the PMM as an nice overview
#'
#' @noRd
#' @export
#'
#' @import data.table
#' @import ggplot2
#' @import grid
#'
#' @importFrom RColorBrewer brewer.pal
plot.svd.analysed.pmm.model.data <-
function
(
  x, y,
  ...
)
{
  params <- list(...)
  size <- ifelse(hasArg(size), params$size, 14)
  LDcolors <- RColorBrewer::brewer.pal(unique(x$Virus),"Spectral")
  p1 <-
    ggplot2::ggplot(x, aes(Virus, Readout)) +
    ggplot2::geom_boxplot(
      aes(fill=Virus),
      outlier.shape = NA, outlier.size=NA, notch=F,
      position="identity") +
    ggplot2::scale_y_continuous(limits = quantile(x$Readout, na.rm=T, c(0.1, 0.9))) +
    ggplot2::facet_grid(. ~ Screen ) +
    ggplot2::scale_fill_manual(values=LDcolors) +
    ggplot2::theme_bw()
  p2 <-
    ggplot2::ggplot(x, aes(Readout, ..density..)) +
    ggplot2::geom_histogram(aes(fill=Virus), bins=100) +
    ggplot2::ylab("Density") +
    ggplot2::facet_grid(Screen ~ Virus) +
    ggplot2::scale_fill_manual(values=LDcolors) +
    ggplot2::theme_bw()

  .multiplot(p1, p2)
}

#' Plots the graph of network diffusion using 1-NN
#'
#' @noRd
#' @export
#'
#' @import data.table
#' @import igraph
#' @importFrom graphics plot legend
plot.svd.diffused.pmm <- function(x, y, ...)
{
   obj <- x$graph.info$graph
   graphics::plot.new()
   V(obj)$size = igraph::degree(obj)
   deg <- igraph::degree(obj)
   size <- deg
   size[deg < 3] <- 15
   size[deg >= 3] <- 20
   size[deg > 5] <- 25
   ad <- igraph::get.adjacency(obj)
   ad[ad >= 1] <- 1
   obj <- igraph::graph_from_adjacency_matrix(ad, mode="undirected")
   blue.genes <- V(x$graph.info$graph)[which(V(x$graph.info$graph)$color == "lightblue")]
   orange.genes <- V(x$graph.info$graph)[which(V(x$graph.info$graph)$color == "orange")]
   V(obj)$color[V(obj) %in% blue.genes] <- "lightblue"
   V(obj)$color[V(obj) %in% orange.genes] <- "orange"
   E(obj)$width <- 2
   graphics::plot(obj, vertex.size=size,layout =  layout.kamada.kawai,
                  vertex.label.family="Helvetica",
                  edge.curved=-.01)
   graphics::legend("topright", legend=x$graph.info$legend, col=x$graph.info$colors,
          pch=19, cex=1.05)
}

#' Plots several plots in one
#'
#' @noRd
#' @import grid
.multiplot <-
  function
(
  ...,
  plotlist=NULL,
  cols=2,
  layout=NULL
)
{
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    return(plots[[1]])
  } else {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
    }
  }
}

