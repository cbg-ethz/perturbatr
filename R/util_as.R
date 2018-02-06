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
#' @import data.table
methods::setAs(
  "data.table",
  "perturbation.data",
  function(from)
  {
    return(methods::new("perturbation.raw.data", .data=from))
	}
)

#' @noRd
#' @import methods
#' @import data.table
methods::setAs(
  "data.frame",
  "perturbation.data",
  function(from)
  {
    methods::as(data.table::as.data.table(from), "perturbation.data")
  }
)
