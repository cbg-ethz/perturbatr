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


#' Plot a \code{knockout.lmm.analysed}  object
#'
#' @export
#' @method plot knockout.lmm.analysed
#' @import ggplot2
#' @import data.table
#' @importFrom dplyr filter
#'
#' @param x  the object to be plotted
#' @param size  size of letters
#' @param ...  additional parameters
#'
#' @return returns a list of plots
#'
plot.knockout.lmm.analysed <- function(x, size=10, ...)
{
  pl <- .plot.knockout.lmm.analysed (x@.gene.hits,
                                     main="Gene effects",
                                     size=size, ...)
  pl2 <-
    .plot.knockout.lmm.analysed(x@.gene.pathogen.hits, main="", size=size) +
    ggplot2::facet_wrap(
      ~ Virus,
      ncol=ceiling(length(unique(x@.gene.pathogen.hits$Virus))/2))
  pl3 <- .plot.effect.matrices.knockout.analysed.lmm(x)
  pl4 <- .plot.hit.counts(x)
  pl5 <- .plot.vulcano(x)
  return(list(gene.effect.barplot=pl,
              gene.pathogen.effect.barplot=pl2,
              gene.pathogen.effect.matrix=pl3,
              gene.pathogen.hit.counts=pl4,
              gene.qval.volcano.plot=pl5))
}

#' @noRd
#' @import data.table
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom methods hasArg
.plot.knockout.lmm.analysed  <- function(x, main, size, ...)
{
  pars <- list(...)
  if ("Virus" %in% colnames(x))
  {
    x <- dplyr::filter(x, Control == 0) %>%
      .[order(abs(Effect), decreasing=TRUE), .SD[1:25], by=Virus] %>%
      dplyr::filter(!is.na(GeneSymbol))
  }
  else
  {
    x <- x[order(abs(Effect), decreasing=TRUE), .SD[1:25]] %>%
      dplyr::filter(!is.na(GeneSymbol), !is.na(Effect))
  }
  x.pos.range <- max(abs(x$Effect))
  x.lim  <- c(-x.pos.range, x.pos.range) + c(-x.pos.range, x.pos.range)/5
  LDcolors <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
  pl <-
    ggplot2::ggplot(x) +
    ggplot2::geom_bar(ggplot2::aes(x=GeneSymbol, y=abs(Effect), fill=Effect),
                      stat="identity") +
    ggplot2::scale_fill_gradient2(low=LDcolors[1],
                                  high=LDcolors[11],
                                  na.value=LDcolors[6],
                                  name="Effect") +
    ggplot2::ylab("Gene effect strength") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y=element_blank(),
                   axis.ticks=element_blank(),
                   text = element_text(size = size, family = "Helvetica"),
                   axis.text.x = element_text(angle=45,
                                              size = size,
                                              family = "Helvetica"),
                   strip.text=element_text(face=x$font))+
    ggplot2::coord_polar() +
    ggplot2::ggtitle(main)

  pl
}

#' @noRd
#' @import data.table
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tidyr gather
.plot.effect.matrices.knockout.analysed.lmm <- function(x, ...)
{

  effect.matrices <- .effect.matrices(x)
  ge <- effect.matrices$gene.effects %>%
    .[order(-abs(Effect))]  %>%
    .[1:25]
  gpe <-  effect.matrices$gene.pathogen.effects %>%
    dplyr::filter(GeneSymbol %in% ge$GeneSymbol) %>%
    tidyr::gather(GeneSymbol)
  colnames(gpe) <- c("GeneSymbol", "Pathogen", "Effect")
  gpe$GeneSymbol <- factor(gpe$GeneSymbol, levels=rev(unique(gpe$GeneSymbol)))
  LDcolors <- rev(RColorBrewer::brewer.pal(11, "Spectral"))

  pl <-
    ggplot2::ggplot(gpe, ggplot2::aes(GeneSymbol, Pathogen)) +
    ggplot2::geom_tile(ggplot2::aes(fill = Effect), colour=LDcolors[1]) +
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::scale_y_discrete(expand = c(0,0)) +
    ggplot2::scale_fill_gradient2(low      = LDcolors[1],
                                  high     = LDcolors[11],
                                  na.value = LDcolors[6],
                                  name     = "Gene virus\neffect") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::theme(text=ggplot2::element_text(size = 8, family = "Helvetica"),
                   aspect.ratio = 2,
                   axis.text.x  = ggplot2::element_text(angle=45,
                                                        hjust = 1,
                                                        size=9),
                   axis.text.y  = ggplot2::element_text(size=9),
                   axis.title   = ggplot2::element_blank(),
                   axis.ticks   = ggplot2::element_blank())
  pl
}

#' @noRd
#' @importFrom dplyr group_by summarize mutate
.plot.hit.counts <- function(x)
{
  obj <- x@.gene.pathogen.effects
  single.res <- obj %>%
    dplyr::mutate(Sign=sign(Effect)) %>%
    dplyr::group_by(Virus, Sign) %>%
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


#' @noRd
#' @import ggplot2
.plot.vulcano <- function(obj, ...)
{
  gpes <- obj@.gene.effects
  gpes$Qval[is.na(gpes$Qval)] <- 1
  x      <- gpes$Effect
  y      <- -log10(gpes$Qval + 0.0001)
  ctrls  <- gpes$Control
  genes  <- gpes$GeneSymbol
  effect.thresh <- obj@.params$effect.size
  sig.thresh    <- -log10(obj@.params$qval.threshold)

  colors <- rep("grey", length(ctrls))
  colors[abs(x) >= effect.thresh & y < sig.thresh] <- "blue"
  colors[which(ctrls == -1)] <- "red"
  colors[which(ctrls == 1)]  <- "green"

  xlim         <- range(x[is.finite(x)])
  ylim         <- range(y[is.finite(y)])
  hits.idx     <- abs(x) >= effect.thresh & y < sig.thresh
  pos.ctrl.idx <- which(ctrls == 1)
  neg.ctrl.idx <- which(ctrls == -1)
  no.hit.idx   <- which(colors == "grey")
  nohits       <- data.frame(Effect=x[no.hit.idx],
                             sig=y[no.hit.idx], Gene=genes[no.hit.idx])
  hits         <- data.frame(Effect=x[hits.idx],
                             sig=y[hits.idx], Gene=genes[hits.idx])
  pos.ctrls    <- data.frame(Effect=x[pos.ctrl.idx],
                             sig=y[pos.ctrl.idx], Gene=genes[pos.ctrl.idx])
  neg.ctrls    <- data.frame(Effect=x[neg.ctrl.idx],
                             sig=y[neg.ctrl.idx], Gene=genes[neg.ctrl.idx])
  pl <-
    ggplot2::ggplot() +
    ggplot2::xlim(xlim) +
    ggplot2::ylim(ylim) +
    ggplot2::xlab("Readout") +
    ggplot2::ylab("-log(q-value)")
  if (length(nohits$Effect) != 0)
    pl <- pl +
    ggplot2::geom_point(aes(x=nohits$Effect, y=nohits$sig),
                        col="grey",size=.5)
  if (length(hits$Effect) != 0)
    pl <- pl +
    ggplot2::geom_point(aes(x=hits$Effect, y=hits$sig, color="Hit"),
                        size=1.5, alpha=.75)
  if (sig.thresh > 0)
    pl <- pl +
    ggplot2::geom_hline(yintercept=sig.thresh, alpha=.75, linetype="dashed") +
    ggplot2::geom_text(
      aes(xlim[1], sig.thresh,
          label = paste("Significance threshold:", sig.thresh),
          vjust = -.25, hjust=0), size = 4)
  if (effect.thresh > 0)
  {
    pl <- pl +
      ggplot2::geom_vline(xintercept=effect.thresh,
                          alpha=.75, linetype="dotdash") +
      ggplot2::geom_vline(xintercept=-effect.thresh,
                          alpha=.75, linetype="dotdash") +
      ggplot2::geom_text(aes(effect.thresh, ylim[2],
                             label = paste("Effect threshold:", effect.thresh),
                             vjust = 0, hjust=-.01), size = 4)
  }
  if (length(pos.ctrls$Effect) != 0)
  {
    pl <- pl +
      ggplot2::geom_point(aes(x=pos.ctrls$Effect,
                              y=pos.ctrls$sig,
                              col="Positive control"), size=1.5, alpha=.75)
  }
  if (length(neg.ctrls$Effect) != 0)
  {
    pl <- pl +
      ggplot2::geom_point(aes(x=neg.ctrls$Effect,
                              y=neg.ctrls$sig,
                              col="Negative control"), size=1.5, alpha=.75)
  }
  if (length(pos.ctrls$Effect) != 0)
  {
    pl <- pl +
      ggplot2::geom_text(aes(x=pos.ctrls$Effect,
                             y=pos.ctrls$sig,
                             label=pos.ctrls$Gene),
                         hjust=0, vjust=0,
                         check_overlap=TRUE,
                         nudge_x=0.05, size=4)
  }
  if (length(neg.ctrls$Effect) != 0)
  {
    pl <- pl +
      ggplot2::geom_text(aes(x=neg.ctrls$Effect,
                             y=(neg.ctrls$sig),
                             label=neg.ctrls$Gene),
                         hjust=0, vjust=0,
                         check_overlap=TRUE,
                         nudge_x=0.05, size=4)
  }
  pl <- pl +
    ggplot2::scale_color_manual(name="Points", values=c("blue", "red", "green")) +
    ggplot2::guides(color=guide_legend(title=NULL)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text=element_text(size=12),
                   axis.title=element_text(size=14,face="bold"),
                   legend.text=element_text(size=14))
  pl
}


#' Plot a \code{knockout.hyper.analysed} object
#'
#' @export
#' @import data.table
#' @method plot knockout.hyper.analysed
#'
#' @param x  the object to plot
#' @param size  size of the text
#' @param ...  additional parameters
#'
#' @return returns a plot object
plot.knockout.hyper.analysed <- function(x, size=10, ...)
{
  .plot.knockout.analysed(x, ...)
}

#' Plot a \code{knockout.tstatistic.analysed} object
#'
#' @export
#' @import data.table
#' @method plot knockout.tstatistic.analysed
#'
#' @param x  the object to plot
#' @param size  size of the text
#' @param ...  additional parameters
#'
#' @return returns a plot object
plot.knockout.tstatistic.analysed <- function(x, size=10, ...)
{
  .plot.knockout.analysed(x, size, ...)
}

#' @noRd
#' @import data.table
#' @import ggplot2
#' @importFrom dplyr group_by summarize mutate filter
.plot.knockout.analysed <- function(x, size, ...)
{
  df <- x@.gene.hits[order(abs(MeanEffect), decreasing=TRUE), .SD[1:25]] %>%
    dplyr::filter(!is.na(GeneSymbol), !is.na(MeanEffect))
  x.pos.range <- max(abs(x@.gene.hits$MeanEffect))
  x.lim  <- c(-x.pos.range, x.pos.range) + c(-x.pos.range, x.pos.range)/5

  LDcolors <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
  pl <- ggplot2::ggplot(df) +
    ggplot2::geom_bar(ggplot2::aes(x=GeneSymbol,y=abs(MeanEffect),
                                   fill=MeanEffect),
                      stat="identity") +
    ggplot2::scale_fill_gradient2(low=LDcolors[1], high=LDcolors[11],
                                  na.value=LDcolors[6],
                                  name="Effect") +
    ylab("Effect") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y=element_blank(),
                   axis.ticks=element_blank(),
                   text = element_text(size = size, family = "Helvetica"),
                   axis.text.x = element_text(angle=45, size = size,
                                              family = "Helvetica"),
                   strip.text = element_text()) +
    ggplot2::coord_polar()

  pl
}
