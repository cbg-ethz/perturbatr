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

#' @export
#' @import ggplot2
#' @import data.table
#' @importFrom RColorBrewer brewer.pal
#' @importFrom methods hasArg
#' @method plot svd.concordance
plot.svd.concordance <- function(x, y, ...)
{
  params <- list(...)
  type <- ifelse(methods::hasArg(type), params$type, "jaccard")
  size <- ifelse(methods::hasArg(size), params$size, 14)
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
    ggplot2::ggtitle(ifelse(methods::hasArg("main"), params$main, "")) +
    ggplot2::guides(fill = guide_colorbar(
      title.position = "top", title.hjust = 0.5))

  pl
}

#' @export
#' @import data.table
#' @import ggplot2
#' @importFrom stats cor
#' @importFrom methods hasArg
#' @method plot svd.replicates
plot.svd.replicates <- function(x, y, ...)
{
  params <- list(...)
  meth <- ifelse(methods::hasArg(method), params$method, "raw")
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

#' @export
#' @import data.table
#' @import ggplot2
#' @method plot svd.replicate
plot.svd.replicate <- function(x, y, ...)
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

#' Plot a \code{knockout.plate}
#'
#'   @export
#' @import data.table
#' @import ggplot2
#' @importFrom dplyr full_join
#' @importFrom methods hasArg
#' @method plot knockout.plate
plot.knockout.plate <- function(x, y,
                                show.controls=FALSE,
                                show.gene.names=FALSE,
                                ylab="Row idx", xlab="Column idx",
                                main="", ...)
{

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
  if (show.controls)
  {
    ctrl <- df$Control
    lwd <- ctrl
    lwd[lwd != 0] <- 1
    pl <- pl + ggplot2::geom_tile(aes(fill=Readout), color="black", lwd=lwd)
  }
  else
  {
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
    ggplot2::ggtitle(main) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size = 12, family = "Helvetica"),
                   aspect.ratio=.75)
  if (show.gene.names) pl <- pl + ggplot2::geom_text(aes(label=df$Gene), size=3)
  pl
}

#' @export
#' @import data.table
#' @import ggplot2
#' @importFrom dplyr select mutate group_indices filter summarize group_by
#' @importFrom tidyr gather
#' @method plot svd.quality
plot.svd.quality <- function(x, y, ...)
{
  # plot the raw plate values as boxplot
  qual <- x$data
  grps <- dplyr::group_indices(qual, Virus, Screen, Library,
                               ScreenType, ReadoutType,
                               Replicate, Plate)
  df   <- dplyr::mutate(qual, Plate=grps) %>%
    dplyr::select(Virus, Screen, Readout, Plate) %>%
    dplyr::mutate(Plate=as.factor(Plate)) %>%
    dplyr::filter(!is.na(Readout), !is.infinite(Readout))
  pl <-
    ggplot2::ggplot(df, ggplot2::aes(Plate, Readout)) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_grid(Virus ~ Screen) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x=ggplot2::element_blank()) +
    ggplot2::ggtitle("Plate readouts")
  # plot control densities
  df  <- dplyr::select(qual, Virus, Screen, Readout, Control) %>%
    dplyr::filter(!is.na(Readout),
                  !is.infinite(Readout),
                  !is.na(Control),
                  Control != 0) %>%
    dplyr::mutate(Control = as.character(Control))
  data.table::setDT(df)[Control == "-1", Control := "Negative control"]
  data.table::setDT(df)[Control == "1", Control := "Positive control"]
  data.table::setDT(df)[Control == "0", Control := "Normal"]
  pl2 <- ggplot2::ggplot() + ggplot2::theme_bw()
  if (nrow(df) != 0) {
    pl2 <-
      ggplot2::ggplot(df) +
      ggplot2::geom_histogram(
        ggplot2::aes(x=Readout, y=..density.., color=Control), alpha=.5,
        position="identity", bins=42) +
      ggplot2::facet_grid(Virus ~ Screen) +
      ggplot2::theme_bw() +
      ggplot2::ylab("Density") +
      ggplot2::ggtitle("Control densities")
  }
  # plot z factor and ssmd per plate as boxplot
  qual <- x$quality$plate.quality %>% ungroup
  df <- dplyr::select(qual, Virus, Screen, z.fac.control, ssmd) %>%
    tidyr::gather(key, value, z.fac.control, ssmd) %>%
    dplyr::filter(!is.na(value), !is.infinite(value))
  pl3 <- ggplot2::ggplot() + ggplot2::theme_bw()
  if (nrow(df) != 0) {
    pl3 <-
      ggplot2::ggplot(df) +
      ggplot2::geom_boxplot(ggplot2::aes(key, value), outlier.shape=NA) +
      ggplot2::scale_y_continuous(limits = quantile(df$value, c(0.1, 0.9))) +
      ggplot2::facet_grid(Virus ~ Screen) +
      ggplot2::theme_bw() +
      ggplot2::xlab("Quality metric") +
      ggplot2::ylab("Readout") +
      ggplot2::scale_x_discrete(breaks=c("ssmd", "z.fac.control"),
                                labels=c("SSMD", "Z-factor")) +
      ggplot2::ggtitle("Plate quality measures")
  }
  # plot positive control and negative control values
  qual <- x$data
  df   <- dplyr::mutate(qual, Plate=grps) %>%
    dplyr::select(Readout, Virus, Screen, Plate, Control) %>%
    dplyr::filter(!is.na(Control), Control != 0) %>%
    dplyr::group_by(Virus, Screen, Plate, Control) %>%
    dplyr::summarize(Readout=mean(Readout, na.rm=T)) %>% ungroup %>%
    dplyr::mutate(Plate=as.factor(Plate), Control=as.factor(Control))
  pl4 <- ggplot2::ggplot() + ggplot2::theme_bw()
  if (nrow(df) != 0) {
    pl4 <-
      ggplot2::ggplot(df) +
      ggplot2::geom_point(ggplot2::aes(x=Plate, y=Readout, color=Control)) +
      ggplot2::theme_bw() +
      ggplot2::geom_line(ggplot2::aes(x=Plate, y=Readout)) +
      ggplot2::facet_grid(Virus ~ Screen) +
      ggplot2::theme(axis.text.x=ggplot2::element_blank()) +
      ggplot2::scale_color_discrete(breaks=c("-1", "1"),
                                    labels=c("Negative control",
                                             "Positive control")) +
      ggplot2::ggtitle("Plate controls")
  }
  .multiplot(plotlist=list(pl, pl2, pl3, pl4),cols=2)
}
