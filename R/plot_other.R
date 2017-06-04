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


#' Plot two \code{knockout.replicate} objects
#'
#' @description Scatter two \code{knockout.replicate} objects against each other
#'
#' @export
#' @method plot knockout.replicate
#' @import data.table
#' @import ggplot2
#' @param x  a \code{knockout.replicate} object
#' @param y  another \code{knockout.replicate} object
#' @param method either \code{Scatterplot} or \code{QQ-plot}
plot.knockout.replicate <- function(x, y, method=c("Scatterplot", "QQ-plot"))
{
  method <- match.arg(method)
  if (missing(x) | missing(y)) stop("Please provide two arguments (x and y)!")

  x@.data <- x@.data %>% .[order(Replicate, Plate, RowIdx, ColIdx)]
  y@.data <- y@.data %>% .[order(Replicate, Plate, RowIdx, ColIdx)]

  df <- data.frame(X=x@.data$Readout, Y=y@.data$Readout)
  if (method == "QQ-plot")
  {
    df$X <- sort(df$X)
    df$Y <- sort(df$Y)
  }

  idx <- which(!is.na(df$X) & !is.na(df$Y))
  corr <- format(cor(df$X[idx], df$Y[idx], method="spearman"), digits=2, width=1)
  pl  <- ggplot2::ggplot(df) +
         ggplot2::geom_point(ggplot2::aes(x=X, y=Y), pch=16) +
         ggplot2::xlab("Replicate 1") +
         ggplot2::ylab("Replicate 2") +
         ggplot2::theme_bw() +
        ggplot2::theme(text = ggplot2::element_text(size=20), aspect.ratio=.75) +
        ggplot2::ggtitle(bquote(paste(.(method), " (", rho == .(corr), ")")))

  pl
}

#' Plot a \code{knockout.plate} object
#'
#' @description Plot a \code{knockout.plate} object on a 2D grid
#'
#' @export
#' @import data.table
#' @import ggplot2
#' @importFrom dplyr full_join
#' @importFrom methods hasArg
#' @method plot knockout.plate
#'
#' @param x  the object to plot
#' @param show.controls  show which wells are controls
#' @param show.gene.names  show the gene names for every well
#' @param xlab  a title for the xlab
#' @param ylab  a title for the ylab
#' @param main  the title of the plot
#' @param axis.text.size  text size of axis labels
#' @param gene.text.size  the size of the gene names within each well
#' @param ...  additional params
plot.knockout.plate <- function(x,
                                show.controls=FALSE,
                                show.gene.names=FALSE,
                                ylab="Row idx",
                                xlab="Column idx",
                                main="",
                                axis.text.size=12,
                                gene.text.size=3,
                                ...)
{
  mat <- plate.matrix(x)
  mat$genes[is.na(mat$genes)] <- ""
  dr  <- data.table::melt(mat$readout)
  di  <- data.table::melt(mat$idx)
  dg  <- data.table::melt(mat$genes)
  df  <- dplyr::full_join(dr, di, by=c("Var1", "Var2"))
  df  <- dplyr::full_join(df, dg, by=c("Var1", "Var2"))
  colnames(df) <- c("Row", "Column", "Readout", "Control", "Gene")
  pl <-
    ggplot2::ggplot(df, ggplot2::aes(x=Column, y=rev(Row)))
  if (show.controls)
  {
    ctrl <- df$Control
    lwd <- ctrl
    lwd[lwd != 0] <- 1
    pl <- pl + ggplot2::geom_tile(ggplot2::aes(fill=Readout),
                                  color="black",
                                  lwd=lwd)
  }
  else
  {
    pl <- pl + ggplot2::geom_tile(ggplot2::aes(fill=Readout),
                                  color="black")
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
    ggplot2::theme(text = ggplot2::element_text(
      size = axis.text.size, family = "Helvetica"), aspect.ratio=.75)
  if (show.gene.names)
    pl <- pl + ggplot2::geom_text(ggplot2::aes(label=df$Gene),
                                  size=gene.text.size)
  pl
}

#' Plot a \code{knockout.quality} object
#'
#' @description Plot a \code{knockout.quality} object on a 2D grid
#'
#' @export
#' @import data.table
#' @import ggplot2
#' @importFrom dplyr select mutate group_indices filter summarize group_by
#' @importFrom tidyr gather
#' @method plot knockout.quality
plot.knockout.quality <- function(x, axis.text.size=12)
{
  # plot the raw plate values as boxplot
  qual <- x@.data
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
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   text = ggplot2::element_text(
                     size = axis.text.size, family = "Helvetica"), aspect.ratio=.75) +
    ggplot2::ggtitle("Plate readouts")

  # plot control densities
  df  <- dplyr::select(qual, Virus, Screen, Readout, Control) %>%
    dplyr::filter(!is.na(Readout),
                  !is.infinite(Readout),
                  !is.na(Control),
                  Control != 0) %>%
    dplyr::mutate(Control = as.character(Control))
  data.table::setDT(df)[Control == "-1", Control := "negative"]
  data.table::setDT(df)[Control == "1",  Control := "positive"]
  data.table::setDT(df)[Control == "0",  Control := "Normal"]
  pl2 <- ggplot2::ggplot() + ggplot2::theme_bw()
  if (nrow(df) != 0)
  {
    pl2 <-
      ggplot2::ggplot(df) +
      ggplot2::geom_histogram(
        ggplot2::aes(x=Readout, y=..density.., color=Control), alpha=.5,
        position="identity", bins=42) +
      ggplot2::facet_grid(Virus ~ Screen) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     text = ggplot2::element_text(
                       size = axis.text.size, family = "Helvetica"), aspect.ratio=.75) +
      ggplot2::ylab("Density") +
      ggplot2::ggtitle("Control densities")
  }

  # plot z factor and ssmd per plate as boxplot
  qual <- x@.quality$plate.quality %>% ungroup
  df <- dplyr::select(qual, Virus, Screen, z.fac.control, ssmd) %>%
    tidyr::gather(key, value, z.fac.control, ssmd) %>%
    dplyr::filter(!is.na(value), !is.infinite(value))
  pl3 <- ggplot2::ggplot() + ggplot2::theme_bw()

  if (nrow(df) != 0)
  {
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
  qual <- x@.data
  df   <- dplyr::mutate(qual, Plate=grps) %>%
    dplyr::select(Readout, Virus, Screen, Plate, Control) %>%
    dplyr::filter(!is.na(Control), Control != 0) %>%
    dplyr::group_by(Virus, Screen, Plate, Control) %>%
    dplyr::summarize(Readout = mean(Readout, na.rm=T)) %>% ungroup %>%
    dplyr::mutate(Plate      = as.factor(Plate), Control=as.factor(Control))
  pl4 <- ggplot2::ggplot() + ggplot2::theme_bw()
  if (nrow(df) != 0)
  {
    pl4 <-
      ggplot2::ggplot(df) +
      ggplot2::geom_point(ggplot2::aes(x=Plate, y=Readout, color=Control)) +
      ggplot2::theme_bw() +
      ggplot2::geom_line(ggplot2::aes(x=Plate, y=Readout)) +
      ggplot2::facet_grid(Virus ~ Screen) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     text = ggplot2::element_text(
                       size = axis.text.size, family = "Helvetica"), aspect.ratio=.75) +
      ggplot2::scale_color_discrete(breaks = c("-1", "1"),
                                    labels = c("negative",
                                               "positive")) +
      ggplot2::ggtitle("Plate controls")
  }

  list(pl, pl2, pl3, pl4)
}
