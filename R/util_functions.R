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
.remove.suffix <- function(s)
{
  if (length(grep("\\.[[:alpha:]]+$", s, perl=T)) != 0)
  {
    s <- sub("\\.[[:alpha:]]+$", "", s)
  }
  s
}

#' @noRd
.is.whole <- function(a)
{
  is.numeric(a) && floor(a) == a
}

#' @noRd
.leuniq <- function(obj)
{
  length(unique(obj))
}

