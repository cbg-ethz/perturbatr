# perturbR: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2018 Simon Dirmeier
#
# This file is part of perturbR
#
# perturbR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# perturbR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with perturbR. If not, see <http://www.gnu.org/licenses/>.



.summarize.two.sided <- . %>%
  dplyr::summarize(HitRatio   = (sum(Hit == TRUE, na.rm=TRUE) / n()),
                   Pval       = metap::sumlog(Pval)$p,
  								 Qval       = metap::sumlog(Qval)$p,
  								 MeanEffect = mean(Readout, na.rm=TRUE),
  								 MaxEffect  = max(Readout, na.rm=TRUE),
  								 MinEffect  = min(Readout, na.rm=TRUE),
  								 MinPval    = min(Pval, na.rm=TRUE),
  								 MinQval    = min(Qval, na.rm=TRUE),
  								 AllPval=paste(sprintf("%03f", Pval), collapse=","),
  								 AllQval=paste(sprintf("%03f", Qval), collapse=",")) %>%
	ungroup
