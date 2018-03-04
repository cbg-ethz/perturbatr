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
.leuniq <- function(obj)
{
  length(unique(obj))
}


#' @noRd
.check <- function(object, cols)
{
  ps   <- paste(cols, collapse=", " )
  coln <- colnames(object)
  if (!all(cols %in% coln))
  {
    msg <- paste("Your data needs cols:", ps)
    msg <- paste(msg, "\nYou have:", paste0(coln, collapse=", "))
    stop(msg)
  }
}

#' @noRd
check.normalized <- function(obj)
{
  if (dataType(obj) != .dataTypes()$NORMALIZED)
      stop("Please preprocess your data before analysing it.")
}
