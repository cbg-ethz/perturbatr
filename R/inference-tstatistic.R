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


#' @title Calculate statistics based on the t-test to analyse the data.
#'
#' @description For this you should use a standardization before.
#'
#' @export
#' @docType methods
#' @rdname t_statistic-methods
#'
#' @import data.table
#'
#' @param obj  the data to be analysed
#' @param mu  side to which the mean of a siRNA is compared to
#' @param padjust  multiple testing correction method
#' @param ...   additional params
setGeneric(
  "t.statistic",
  function(obj,
           mu=c(0, "scrambled", "control"),
           padjust=c("BH", "bonferroni"),
           ...)
  {
    standardGeneric("t.statistic")
  },
  package="knockout"
)

#' @rdname t_statistic-methods
#' @aliases t.statistic,knockout.lmm.data-method
#' @import data.table
setMethod(
  "t.statistic",
  signature=signature(list(obj="knockout.data")),
  function(obj,
           mu=c(0, "scrambled", "control"),
           padjust=c("BH", "bonferroni"),
           ...)
  {
    .check.data(obj)
    dat <- obj@.data
    res <- .t.statistic(dat,
                       mu=match.arg(mu),
                       padjust=match.arg(padjust),
                       ...)
    ret <- new("knockout.analysed",
               .inference=.inference.types()$T.TEST,
               .data=res)
    ret
  }
)

#' @noRd
#' @import data.table
.t.statistic <- function(obj, mu, padjust, ...)
{
  message(paste("Correcting with ", padjust, "!", sep=""))
  message(paste("Taking", mu, "for t-test mu!"))

  ret <- dplyr::group_by(obj, Virus, Screen, Library, ReadoutType, ScreenType,
                         Design, Cell) %>%
    dplyr::mutate(grp=.GRP) %>%
    ungroup
  grps <- unique(ret$grp)
  res <- do.call(
    "rbind",
    lapply(
      grps,
      function (g)
      {
        grp.dat <- dplyr::filter(ret, grp==g) %>% ungroup
        message(paste("Doing grp: ", paste(grp.dat$Virus[1],
                                           grp.dat$Screen[1],
                                           grp.dat$ScreenType[1],
                                           grp.dat$ReadoutType[1],
                                           grp.dat$Cell[1],
                                           grp.dat$Design[1],
                                           grp.dat$Library[1],
                                           sep=", ")))
        fr <- .do.t.statistic(grp.dat, padjust, mu)
        fr
      }
    )
  )
  res
}


#' @noRd
#' @import data.table
#' @importFrom dplyr mutate group_by summarize ungroup
#' @importFrom stats p.adjust
#' @importFrom asserthat assert_that
.do.t.statistic <- function(obj, padjust, mu)
{
  res <- obj %>% ungroup %>%
    dplyr::group_by(Plate) %>%
    dplyr::mutate(grp=.GRP)
  grps <- unique(res$grp)
  ret  <- do.call(
    "rbind",
    lapply(
      grps,
      function(g)
      {
        grp.dat <- dplyr::filter(res, grp==g) %>% ungroup
        fr      <- .t.statisic.plate(grp.dat, mu)
        fr
      }
    )
  )
  data.table::setDT(ret)[,Pval := as.numeric(Pval)]
  data.table::setDT(ret)[,Qval := p.adjust(Pval, method=padjust)]
  assertthat::assert_that(all(order(ret$Pval) == order(ret$Qval)))
  ret
}

#' @noRd
#' @import data.table
#' @importFrom dplyr filter group_by mutate summarize
.t.statisic.plate <- function(obj, mu)
{
  ret <- obj %>% ungroup
  mu <- .set.mu(ret, mu)
  ret <- dplyr::group_by(ret, Virus, Screen, Library,
                         ReadoutType, ScreenType, Cell, Design,
                         GeneSymbol, Entrez, Plate, Control,
                         RowIdx, ColIdx, siRNAIDs) %>%
    dplyr::summarize(Pval=.t.test(GeneSymbol, Readout, mu),
                     Readout=mean(Readout, na.rm=T))  %>%
    ungroup
  ret
}

#' @noRd
#' @import data.table
#' @importFrom dplyr filter
.set.mu <- function(ret, mu.gene)
{
  if (mu.gene != "0")
  {
    cont.tab <- dplyr::filter(ret, Control  == -1)
    if (tolower(mu.gene) == "scrambled")
      cont.tab <- dplyr::filter(cont.tab, tolower(GeneSymbol) == "scrambled")
    if (nrow(cont.tab) == 0) stop("No controls found for criteria!")
    if (nrow(cont.tab) < 3)
      stop("Less than three controls found better take mu=0!")
    mu <- cont.tab$Readout
  }
  else if (mu.gene == "0")
  {
    mu <- 0
  }
  else message(paste("\t..taking mu=", mu, ".",sep=""))
  if (!is.numeric(mu)) stop("Please provide a numeric mu!")

  mu
}

#' @noRd
#' @importFrom stats t.test
.t.test <- function(g, val, mu)
{
  tst <- 1
  if (length(val) < 3)
  {
    warning(paste("<3 values provided for" , g, " -> returning 1"))
  }
  else
  {
    tryCatch ({
      if (mu == 0) {
        tst <- stats::t.test(val, mu=0, alternative="two.sided")$p.value
      }
      else {
        tst <- stats::t.test(val, y=mu, alternative="two.sided")$p.value
      }

    }, warning = function(war)
      { warning(paste(war, " -> setting one for", g)); },
       error = function(err)
      { warning(paste(err, " -> setting one for", g)); })
  }
  tst
}

