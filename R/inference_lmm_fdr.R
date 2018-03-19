# perturbatr: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2018 Simon Dirmeier
#
# This file is part of perturbatr
#
# perturbatr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# perturbatr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with perturbatr. If not, see <http://www.gnu.org/licenses/>.


#' @noRd
#' @import tibble
#' @importFrom dplyr mutate select left_join
#' @importFrom tidyr spread
#' @importFrom lme4 ranef
#' @importFrom rlang .data
ge.fdrs <- function(md, ref, bootstrap.cnt, frm)
{
  if (!is.numeric(bootstrap.cnt) | bootstrap.cnt < 10)
  {
    dt <- tibble::tibble(
      GeneSymbol=ref$gene.effects$GeneSymbol, Qval=NA_real_)
    return(list(ret=dt, btst=FALSE))
  }

  message("Bootstrapping ... this might take a while.")
  li <- list()
  i  <- ctr <- 1
  mistrial.cnt <- bootstrap.cnt * 10
  while (i <= bootstrap.cnt)
  {
    ctr <- ctr + 1
    tryCatch({
        bt.sample <- bootstrap(md, .data$Condition, .data$Perturbation)
        lmm.fit   <- .hm.fit(bt.sample, frm)
        re        <- .ranef(lmm.fit)
        da <- tibble::tibble(
            bootstrap  = paste0("Bootstrap_", sprintf("%03i", i)),
            Effect     = re$gene.effects$Effect,
            GeneSymbol = re$gene.effects$GeneSymbol)
        li[[i]] <- da
        i <- i + 1
    }, error = function(e) stop(e),
      warning = function(e) {
      warning(e)
    })
    if (ctr == mistrial.cnt)
        stop(paste0("Breaking bootstrapping after ",
                    mistrial.cnt ," mis-trials!"))
  }

  btst.dat <- dplyr::bind_rows(li)
  fdrs     <- .ge.fdrs(btst.dat, bootstrap.cnt)

  # join bootstrap table with fdr table
  ret  <- dplyr::left_join(
      fdrs, tidyr::spread(btst.dat, .data$bootstrap,
                          .data$Effect), by="GeneSymbol")

  list(ret=ret, btst=TRUE)
}


#' @noRd
#' @import dplyr
#' @importFrom assertthat assert_that
#' @importFrom stats p.adjust
#' @importFrom rlang .data
.ge.fdrs <- function(btst.dat, cnt)
{
  res <- dplyr::group_by(btst.dat, .data$GeneSymbol)
  res <- ungroup(dplyr::do(res, conf.int(.data$Effect, cnt)))
  res <- dplyr::arrange(res, .data$Pval)
  res <- dplyr::mutate(res, "Qval" = p.adjust(.data$Pval, method="BH"))

  res
}


#' @noRd
#' @import tibble
#' @importFrom stats t.test
conf.int <- function(eff, cnt)
{
  t <- stats::t.test(eff, mu=0, na.rm=TRUE)
  tibble::tibble(Mean  = mean(eff, na.rm=TRUE),
                 Pval  = t$p.value,
                 Lower = t$conf.int[1],
                 Upper = t$conf.int[2])
}
