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


#' @aliases rbind,knockout.data-method
#' @import data.table
setMethod(
  "rbind",
  "knockout.data",
  function(...)
  {
    args  <- list(...)
    if (length(args) < 2) return(args[[1]])
    types <- unlist(lapply(args, function(e) e@.type))
    if(all(types == types[1])) stop("Data-types do not agree")
    dat   <- data.table::rbindlist(lapply(args, function(e) e@.data))
    new("knockout.data", .data=dat, .type=types[1])
  }
)
