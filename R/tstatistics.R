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

#' Calculate statistics based on the t-test to analyse the data.
#'
#' For this you should use a standardization before.
#'
#' @export
#' @import data.table
#'
#' @param obj  the data to be analysed
#' @param mu  side to which the mean of a siRNA is compared to
#' @param padjust  multiple testing correction method
#' @param ...   additional params
tstatistic <- function(obj, mu=c(NA, "Scrambled", "control"),
                       padjust=c("BH", "bonferroni"), ...)
{
  UseMethod("tstatistic", obj)
}

#' @noRd
#' @export
#' @import data.table
tstatistic.svd.data <- function(obj, mu=c(NA, "Scrambled", "control"),
                                padjust=c("BH", "bonferroni"), ...)
{
  mu <- match.arg(mu)
  padjust <- match.arg(padjust)
  ret <- .tstatistic(obj, mu=mu, padjust=padjust, ...)
  class(ret) <- c("svd.analysed.tt", "svd.analysed", class(ret))
  invisible(ret)
}

#' @noRd
#' @import data.table
.tstatistic <- function(obj, mu, padjust, ...)
{
  message(paste("Correcting with ", padjust, "!", sep=""))
  message(paste("Taking", mu, "for t-test!"))
  ret <- dplyr::group_by(obj, Virus, Screen, Library,
                         ReadoutType, InfectionType,
                         Design, Cell) %>%
    dplyr::mutate(grp=.GRP) %>%
    ungroup
  grps <- unique(ret$grp)
  res <- do.call(
    "rbind",
    lapply
    (
      grps,
      function (g)
      {
        grp.dat <- dplyr::filter(ret, grp==g) %>% ungroup
        message(paste("Doing grp: ", paste(grp.dat$Virus[1],
                                         grp.dat$Screen[1],
                                         grp.dat$InfectionType[1],
                                         grp.dat$ReadoutType[1],
                                         grp.dat$Cell[1],
                                         grp.dat$Design[1],
                                         grp.dat$Library[1],
                                         sep=", ")))
        fr <- .do.tstatistic(grp.dat, padjust, mu)
        fr
      }
    )
  )
  invisible(res)
}


#' @noRd
#' @import data.table
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr ungroup
#' @importFrom stats p.adjust
.do.tstatistic <- function(obj, padjust, mu)
{
  res <- obj %>% ungroup
  # unfortunately it is needed (but works bcs siRNAs are always on the same plates)
  res <- dplyr::group_by(res, Plate) %>%
    dplyr::mutate(grp=.GRP)
  grps <- unique(res$grp)
  ret <- do.call(
    "rbind",
    lapply
    (
      grps,
      function(g)
      {
        grp.dat <- dplyr::filter(res, grp==g) %>% ungroup
        fr <- .tstatisic.plate(grp.dat, mu)
        fr
      }
    )
  )
  data.table::setDT(ret)[,Pval := as.numeric(Pval)]
  ret <- ret[order(Pval)]
  ret <- dplyr::mutate(ret, Pvalcorr=p.adjust(Pval, method=padjust))
  invisible(ret)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr filter group_by mutate summarize
.tstatisic.plate <- function(obj, mu)
{
  ret <- obj %>% ungroup
  mu <- .set.mu(ret, mu)
  ret <- dplyr::group_by(ret, Virus, Screen, Library,
                         ReadoutType, InfectionType,
                         Cell, Design,
                         GeneSymbol, Entrez, Plate, Control,
                         RowIdx, ColIdx, siRNAIDs) %>%
    dplyr::summarize(Pval=.ttest(GeneSymbol, Readout, mu),
                     Readout=mean(Readout, na.rm=T))  %>%
    ungroup
  ret
}

#' Calculate the mean of the negative controls or the scrambled RNA
#' @noRd
#' @import data.table
#' @importFrom dplyr filter
.set.mu <- function(ret,mu.gene)
{
  if (!is.na(mu.gene))
  {
    cont.tab <- dplyr::filter(ret, Control  == -1)
    if (toupper(mu.gene) == "SCRAMBLED")
      cont.tab <- dplyr::filter(cont.tab, GeneSymbol == "Scrambled")
    if (nrow(cont.tab) == 0) stop("No controls found for criteria!")
    if (nrow(cont.tab) < 3)
      stop(paste("Less than three controls found,",
                 "better standardize data and take mu=0!"))
    mu <- mean(cont.tab$Readout, na.rm=T)
    message(paste("\t..taking mu from all controls(",mu.gene, "): ",
                  mu, sep=""))
  }
  else if (is.na(mu.gene))
  {
    message("\t..taking mu=0.")
    mu <- 0
  }
  else message(paste("\t..taking mu=", mu, ".",sep=""))
  if (!is.numeric(mu)) stop("Please provide a numeric mu!")
  mu
}

#' @noRd
#' @importFrom stats t.test
.ttest <- function(g, val, mu)
{
  tst <- 1
  if (length(val) < 3)
  {
    warning(paste("Only", length(val),
                  "value(s) provided for" , g, "->returning 1"))
  }
  else
  {
    tryCatch ({
      tst <- stats::t.test(val, mu=mu,
                           alternative="two.sided", paired=F)$p.value
    }, warning = function(war)
      { warning(paste(war, " -> setting one for", g)); },
       error = function(err)
      { warning(paste(err, " -> setting one for", g)); }
    )
  }
  invisible(tst)
}

