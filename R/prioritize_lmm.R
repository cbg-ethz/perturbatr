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
#' @importFrom dplyr filter select group_by mutate summarize full_join
.prioritize.lmm <- function(obj,
                            effect.size=0,
                            pval.threshold=0.05,
                            qval.threshold=1)
{
  ge <-
    obj$gene.effects %>%
    dplyr::filter(abs(Effect) >= effect.size)  %>%
    .[order(-abs(Effect))]
  if (!all(is.na(ge$FDR))) ge <- dplyr::filter(ge, FDR <= fdrt)
  gpe <- obj$gene.pathogen.effects %>%
    dplyr::filter(abs(Effect) >= eft)  %>%
    .[order(-abs(Effect))]
  if (!all(is.na(gpe$FDR))) gpe <- dplyr::filter(gpe, FDR <= fdrt)
  list(gene.hits=ge, gene.pathogen.hits=gpe)

}
