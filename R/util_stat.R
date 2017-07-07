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
#' @importFrom stats median
.summarization.method <- function(summ.method)
{
  f <- switch(as.character(summ.method),
              "mean"=base::mean,
              "median"=stats::median,
              "min"=base::min,
              `NA`=NA,
              stop("wrong method given"))
  f
}

#' @noRd
#' @importFrom tibble data_frame
#' @importFrom stats t.test
conf.int <- function(eff, cnt)
{
  # TODO: maybe compute directly on the quantiles for cnt > 1000
  t          <- stats::t.test(eff, mu=0, na.rm=TRUE)
  tibble::data_frame(Mean  = mean(eff, na.rm=TRUE),
                     Pval  = t$p.value,
                     Lower = t$conf.int[1],
                     Upper = t$conf.int[2])
}
