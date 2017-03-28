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
#' @importFrom dplyr filter
#' @importFrom methods hasArg
#' @method plot svd.analysed.pmm
#' @method plot svd.analysed.pmm
plot.svd.analysed.pmm <- function(x, y, ...)
{
  params  <- list(...)
  sig.thresh     <- ifelse(methods::hasArg(sig.thresh),
                           params$fdr.thresh, .2)
  effect.thresh <- ifelse(methods::hasArg(effect.thresh),
                          params$effect.thresh, 0.0)
  log.scale      <- ifelse(methods::hasArg(log.scale),
                           params$log.scale, F)
  gpes <-
    x$gene.pathogen.effects %>%
    dplyr::filter(!is.na(FDR))
  y <- gpes$FDR + 0.0001
  if (log.scale)
  {
    y <- -log(y)
    fdr.thresh <- -log(sig.thresh)
  }
  xl <- ifelse(methods::hasArg(xlab), params$xlab, "Readout")
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

#' @export
#' @import data.table
#' @importFrom dplyr filter
#' @importFrom methods hasArg
#' @method plot svd.analysed
plot.svd.analysed <- function(x, y, ...)
{
  params  <- list(...)
  sig.thresh     <- ifelse(methods::hasArg(sig.thresh), params$sig.thresh, 0.0)
  readout.thresh <- ifelse(methods::hasArg(readout.thresh), params$readout.thresh, 0.0)
  log.scale      <- ifelse(methods::hasArg(log.scale), params$log.scale, F)
  xl <- ifelse(methods::hasArg(xlab), params$xlab, "Readout")
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
.plot.svd.analysed <- function(x, y, ctrls, genes,
                               readout.thresh, sig.thresh,
                               xl, yl, ...)
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
    ggplot2::xlab(xl) +
    ggplot2::ylab(yl)
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
    ggplot2::geom_text(aes(xlim[1], sig.thresh,
                           label = paste("Significance threshold:", sig.thresh),
                           vjust = -.25, hjust=0), size = 4)
  if (readout.thresh > 0)
    pl <- pl +
    ggplot2::geom_vline(xintercept=readout.thresh,  alpha=.75, linetype="dotdash") +
    ggplot2::geom_vline(xintercept=-readout.thresh, alpha=.75, linetype="dotdash") +
    ggplot2::geom_text(aes(readout.thresh, ylim[2],
                           label = paste("Effect threshold:", readout.thresh),vjust = 0, hjust=-.01), size = 4)
  if (length(pos.ctrls$Effect) != 0) pl <- pl +
    ggplot2::geom_point(aes(x=pos.ctrls$Effect,
                            y=pos.ctrls$sig, col="Positive control"), size=1.5, alpha=.75)
  if (length(neg.ctrls$Effect) != 0) pl <- pl +
    ggplot2::geom_point(aes(x=neg.ctrls$Effect,
                            y=neg.ctrls$sig, col="Negative control"), size=1.5, alpha=.75)
  if (length(pos.ctrls$Effect) != 0) pl <- pl +
    ggplot2::geom_text(aes(x=pos.ctrls$Effect,
                           y=pos.ctrls$sig, label=pos.ctrls$Gene), hjust=0, vjust=0, check_overlap=TRUE, nudge_x=0.05, size=4)
  if (length(neg.ctrls$Effect) != 0) pl <- pl +
    ggplot2::geom_text(aes(x=neg.ctrls$Effect,
                           y=(neg.ctrls$sig), label=neg.ctrls$Gene), hjust=0, vjust=0, check_overlap=TRUE, nudge_x=0.05, size=4)
  pl <- pl +
    ggplot2::scale_color_manual(name="Points", values=c("blue", "red", "green")) +
    ggplot2::guides(color=guide_legend(title=NULL)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text=element_text(size=12),
                   axis.title=element_text(size=14,face="bold"),
                   legend.text=element_text(size=14))
  pl
}

#' @export
#' @import data.table
#' @import ggplot2
#' @import grid
#' @importFrom RColorBrewer brewer.pal
#' @importFrom methods hasArg
#' @method plot svd.analysed.pmm.model.data
plot.svd.analysed.pmm.model.data <- function(x, y, ...)
{
  params <- list(...)
  size <- ifelse(methods::hasArg(size), params$size, 14)
  LDcolors <- RColorBrewer::brewer.pal(length(unique(as.character(x$Virus))), "Spectral")
  p1 <-
    ggplot2::ggplot(x, aes(Virus, Readout)) +
    ggplot2::geom_boxplot(
      aes(fill=Virus),
      outlier.shape = NA, outlier.size=NA, notch=F,
      position="identity") +
    ggplot2::scale_y_continuous(limits = quantile(x$Readout, na.rm=T, c(0.01, 0.99))) +
    ggplot2::scale_fill_manual(values=LDcolors) +
    ggplot2::theme_bw()
  p2 <-
    ggplot2::ggplot(x, aes(Readout, ..density..)) +
    ggplot2::geom_histogram(aes(fill=Virus), bins=100) +
    ggplot2::ylab("Density") +
    ggplot2::facet_grid(Virus ~ .) +
    ggplot2::scale_fill_manual(values=LDcolors) +
    ggplot2::theme_bw()

  .multiplot(plotlist=list(p1, p2))
}

#' @export
#' @import ggplot2
#' @method plot svd.analysed.pmm.fdr
plot.svd.analysed.pmm.fdr <- function(x, y, ...)
{
  hits <- x$hist.dat
  # the values used for estimation of the densities (f0 and f1)
  zvals <- data.frame(Z=hits$zvalues)
  # frame of estimated non-null hits (yt),
  # density of mixture(f) and density of nulls (f0)
  densities <- data.frame(X=hits$x, y=hits$yt, f=hits$f, f0=hits$f0)
  pl <-
    ggplot2::ggplot(zvals) +
    ggplot2::geom_histogram(aes(Z), bins = 130, alpha = 0.3) +
    ggplot2::geom_bar(data=densities, aes(x=X, y=y),
                      stat = "identity", color="blue") +
    ggplot2::geom_line(data=densities,
                       aes(x=X, y=f0,  colour="p0*f0"), size=1, linetype=1) +
    ggplot2::geom_line(data=densities,
                       aes(x=X, y=f,  colour="f"), size=1, linetype=2) +
    ggplot2::ylab("Frequency") +
    ggplot2::xlab("Gene-pathogen effect") +
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size=14)) +
    ggplot2::scale_colour_manual(name="Density", values=c("green", "blue"))

  pl
}

