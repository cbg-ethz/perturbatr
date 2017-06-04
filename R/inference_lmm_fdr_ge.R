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
ge.fdrs <- function(md, ref, bootstrap.cnt)
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
                                 Qval=NA_real_)
    return(list(ret=dt, btst=btst))
  }

  li   <- list()
  i    <- 1
  ctr  <- 1
  mistrial.cnt <- bootstrap.cnt * 10
  # TODO: parallelize
  repeat
  {
    ctr <- ctr + 1
    tryCatch(
    {

      bt.sample <- bootstrap(md)
      lmm.fit   <- .lmm(bt.sample)
      re        <- .ranef(lmm.fit)

      da <- data.table::data.table(
        bootstrap  = paste0("Bootstrap_", sprintf("%03i", i)),
        Effect     = re$gene.effects$Effect,
        GeneSymbol = re$gene.effects$GeneSymbol
      )

      li[[i]] <- da
      i <- i + 1

      }, error = function(e) {
        print(paste("Didn't fit:", i, ", error:", e)); i <<- 1000
      }, warning = function(e) {
        print(e)
      }
    )
    if (i > bootstrap.cnt) break
    if (ctr == mistrial.cnt)
      stop(paste0("Breaking after ", mistrial.cnt ," mis-trials!"))
  }
  btst.dat <- data.table::rbindlist(li)
  fdrs     <- .ge.fdrs(btst.dat, bootstrap.cnt)

  # join bootstrap table with fdr table
  ret  <- dplyr::left_join(
    fdrs, tidyr::spread(btst.dat, bootstrap, Effect), by="GeneSymbol")

  list(ret=ret, btst=btst)
}

#' @noRd
#' @importFrom assertthat assert_that
.ge.fdrs <- function(btst.dat, cnt)
{
  res <-
    dplyr::group_by(btst.dat, GeneSymbol) %>%
    dplyr::do(conf.int(.$Effect, cnt)) %>%
    ungroup %>%
    .[order(Pval)] %>%
    dplyr::mutate(Qval=p.adjust(Pval, method="BH"))
  assertthat::assert_that(all(order(res$Pval)  == order(res$Qval)))

  res
}

