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


#' @noRd
#' @import data.table
#' @importFrom dplyr mutate select left_join
#' @importFrom tidyr spread
#' @importFrom lme4 ranef
.ge.fdrs <- function(md, ref, bootstrap.cnt)
{

  if (is.numeric(bootstrap.cnt) & bootstrap.cnt >= 10)
  {
    btst <- TRUE
    message("Bootstrapping...")
  }
  else
  {
    btst <- FALSE
    dt <- data.table::data.table(GeneSymbol=ref$gene.effects$GeneSymbol,
                                 FDR=NA_real_)
    return(list(dt=dt, btst=btst))
  }

  li   <- list()
  i    <- 1
  ctr  <- 1
  mistrial.cnt <- bootstrap.cnt * 10
  repeat
  {
    ctr <- ctr + 1
    tryCatch({
      bt.sample <- bootstrap(md)
      # here use the lmer params
      lmm.fit <- .lmm(bt.sample)
      re      <- .ranef(lmm.fit)
      da <- data.table::data.table(bootstrap=paste0("Bootstrap_", sprintf("%03i", i)),
                                   Effect=re$gene.effects$Effect,
                                   GeneSymbol=re$gene.effects$GeneSymbol)
      li[[i]] <- da
      i <- i + 1
    }, error=function(e) { print(paste("Didn't fit:", i, ", error:", e)); i <<- 1000 },
    warning=function(e) { print(e)})
    if (i > bootstrap.cnt)
      break
    if (ctr == mistrial.cnt)
      stop(paste0("Breaking after ", mistrial.cnt ," mis-trials!"))
  }
  flat.dat <- do.call("rbind", lapply(li, function(e) e))
  # TODO modify that accordingly
  dat <-  flat.dat %>%
    dplyr::group_by(GeneSymbol) %>%
    dplyr::summarise(MeanBootstrap=mean(Effect, na.rm=T),
                     Pval=.ttest(GeneSymbol, Effect, 0)) %>%
    ungroup %>%
    dplyr::mutate(FDR=p.adjust(Pval, method=padj)) %>%
    .[order(FDR)]
  dat <- dplyr::left_join(dat, tidyr::spread(flat.dat, bootstrap, Effect), by="GeneSymbol")
  dat
}
