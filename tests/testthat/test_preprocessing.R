
# knockdown: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2015 - 2016 Simon Dirmeier
#
# This file is part of knockdown
#
# knockdown is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# knockdown is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with knockdown. If not, see <http://www.gnu.org/licenses/>.


context("preprocessing")

data(rnaiscreen)
pl <- plates(rnaiscreen)[1]
pl@.data$PlateIndex <- NULL
plate.data <- methods::as(pl@.data, "knockdown.data")

testthat::test_that("z scoring gives approx standardized gaussian with mean 0", {
  pl.norm <- preprocess(plate.data , normalize="z.score")
  testthat::expect_equal(mean(pl.norm@.data$Readout), 0, tolerance=0.1)
})

testthat::test_that("z scoring gives approx standardized gaussian with var 1", {
  pl.norm <- preprocess(plate.data , normalize="z.score")
  testthat::expect_equal(var(pl.norm@.data$Readout), 1, tolerance=0.1)
})

testthat::test_that("loess normalisation does not change when numcells = NA", {
  loc.plate.data <- plate.data
  loc.plate.data@.data$NumCells <- NA
  pl.norm <- preprocess(loc.plate.data , normalize="loess")
  testthat::expect_equal(pl.norm@.data$Readout,
                         loc.plate.data@.data$Readout,
                         tolerance=0.1)
})

testthat::test_that("bscore normalisation is correct", {
  loc.plate.data <- plate.data
  pl.norm <- preprocess(loc.plate.data , normalize="b.score")
  testthat::expect_equal(
    pl.norm@.data$Readout,
    knockdown:::.bscore.plate(loc.plate.data@.data$RowIdx,
                             loc.plate.data@.data$ColIdx,
                             loc.plate.data@.data$Readout),
    tolerance=0.1)
})
