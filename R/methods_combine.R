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


#' @title Bind multiple perturbation data sets together by row
#'
#' @description Binds multiple \code{PerturbationData} objects together by row.
#'
#' @export
#' @method rbind PerturbationData
#'
#' @import tibble
#' @importFrom dplyr bind_rows
#' @importFrom methods new
#'
#' @param ...  variable number of \code{PerturbationData} objects
#' @return  returns a combined object of class \code{PerturbationData}
#' @examples
#'   data(rnaiscreen)
#'   rbind(rnaiscreen, rnaiscreen)
#'
rbind.PerturbationData <-  function(...)
{
  args  <- list(...)
  if (length(args) < 2) return(args[[1]])
  clazz <- base::unlist(base::lapply(args, function(e) class(e)[1]))
  if(any(clazz != clazz[1])) stop("Data classes do not agree")
  dat   <- dplyr::bind_rows(base::lapply(args, function(e) dataSet(e)))

  methods::new(clazz[1], dataSet=dat)
}
