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
#' @importFrom dplyr group_by summarize ungroup filter select
.prioritize.hyper.statistic <- function(obj,
                                        hit.ratio=0.5,
                                        effect.size=0,
                                        pval.threshold=0.05,
                                        qval.threshold=1)
{
  res <- dplyr::group_by(obj, Virus, Screen, Library,
                         ScreenType, ReadoutType,
                         Design, Cell,
                         GeneSymbol, Entrez) %>%
    dplyr::summarize(HitRatio   = (sum(Hit == TRUE, na.rm=T)/n()),
                     PvalRatio  = (sum(Pval  <= pval.threshold, na.rm=T)/n()),
                     QvalRatio  = (sum(Qval  <= qval.threshold, na.rm=T)/n()),
                     MeanEffect = mean(Readout,na.rm=T),
                     MaxEffect  = max(Readout, na.rm=T),
                     MinEffect  = min(Readout, na.rm=T),
                     MinPval    = min(Pval, na.rm=T),
                     MinQval    = min(Qval, na.rm=T),
                     Pval=paste(sprintf("%.03f", Pval), collapse=","),
                     Qval=paste(sprintf("%.03f", Qval), collapse=",")) %>%
    ungroup %>%
    dplyr::filter(HitRatio >= hit.ratio)
  res
}
