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
#' @import data.table
#' @method plot svd.raw
plot.svd.raw <- function(x, y, ...)
{
  x <- dplyr::filter(x, ReadoutClass=="Readout")
  plot.svd.data(x, ...)
}

#' @export
#' @import ggplot2
#' @import data.table
#' @importFrom dplyr summarize
#' @importFrom dplyr group_by
#' @importFrom tidyr gather
#' @method plot svd.data
plot.svd.data <- function(x, y, ...)
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


#' @export
#' @import data.table
#' @importFrom graphics par
#' @importFrom methods hasArg
#' @method plot svd.plates
plot.svd.plates <- function(x, y, ...)
{
  params <- list(...)
  count <- ifelse(methods::hasArg(count), params$count, NA)
  show.controls <- ifelse(methods::hasArg(show.controls) &
                            is.logical(params$show.controls),
                  params$show.controls, T)
  show.gene.names <- ifelse(methods::hasArg(show.gene.names) &
                              is.logical(params$show.gene.names),
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
                        ", infection:", plate$ScreenType[1],
                        sep="")
    print(plot.svd.plate(plate, main=plate.name,
                         show.controls=show.controls,
                         show.gene.names=show.gene.names))
  }
  graphics::par(ask=F)
}

#' @export
#' @import data.table
#' @import ggplot2
#' @importFrom dplyr full_join
#' @importFrom methods hasArg
#' @method plot svd.plate
plot.svd.plate <- function(x, y, ...)
{
  params <- list(...)
  show.controls <- ifelse(methods::hasArg(show.controls) &
                            is.logical(params$show.controls),
                          params$show.controls, T)
  show.gene.names <- ifelse(methods::hasArg(show.gene.names) &
                              is.logical(params$show.gene.names),
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
    ggplot2::ggtitle(ifelse(methods::hasArg(main), params$main, "")) +
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
    ggplot2::ggplot(df, aes(Plate, Readout)) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_grid(Virus ~ Screen) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x=element_blank()) +
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
                                  labels=c("Negative control",
                                           "Positive control")) +
    ggplot2::ggtitle("Plate controls")

  .multiplot(plotlist=list(pl, pl2, pl3, pl4),cols=2)
}

#' @export
#' @import data.table
#' @import igraph
#' @importFrom graphics plot legend
#' @importFrom methods hasArg
#' @method plot svd.diffused.pmm
plot.svd.diffused.pmm <- function(x, y, ...)
{
   pars <- list(...)
   sz <- ifelse(methods::hasArg(size), pars$size, -1)
   show.labels <- ifelse(methods::hasArg(size), pars$size, -1)
   obj <- x$graph.info$graph
   V(obj)$size = igraph::degree(obj)
   deg <- igraph::degree(obj)
   size <- deg
   size[deg < 3] <- 15
   size[deg >= 3] <- 20
   size[deg > 5] <- 25
   ad <- igraph::get.adjacency(obj)
   ad[ad >= 1] <- 1
   obj <- igraph::graph_from_adjacency_matrix(ad, mode="undirected")
   blue.genes <-
     V(x$graph.info$graph)[which(V(x$graph.info$graph)$color == "lightblue")]
   orange.genes <-
     V(x$graph.info$graph)[which(V(x$graph.info$graph)$color == "orange")]
   V(obj)$color[V(obj) %in% blue.genes] <- "lightblue"
   V(obj)$color[V(obj) %in% orange.genes] <- "orange"
   E(obj)$width <- 2
   graphics::plot.new()
   op <- par(family = "Helvetica", font=2)
   if (sz != -1) size <- rep(sz, length(size))
   graphics::plot(obj, vertex.size=size,layout =  layout.kamada.kawai,
                  vertex.label.family="Helvetica", vertex.label.font=2,
                  edge.curved=-.01)
   graphics::legend("topright",
                    legend=c("Linear mixed model", "Diffusion"),
                    col=x$graph.info$colors,
          pch=19, cex=1.05)
   par(op)
}

#' @noRd
#' @import grid
.multiplot <- function(..., plotlist=NULL, cols=2, layout=NULL)
{
  plots <- c(..., plotlist)
  numPlots <- length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    return(plots[[1]])
  } else {
    grid::grid.newpage()
    grid::pushViewport(
      grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots)
    {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
    }
  }
}
