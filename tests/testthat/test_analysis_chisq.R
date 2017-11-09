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


context("chisq hypothesis tests")

data(rnaiscreen)
rnai.norm <- preprocess(rnaiscreen, normalize="z.score")

testthat::test_that("error is thrown when more than one replicate", {
  testthat::expect_error(
    chisq.statistic(rnai.norm, effect.size=0.01, qval.threshold=1))
})

testthat::test_that("mahalanobis is computed correctly", {
  x <- -5:5
  ord.x <- order(abs(x))
  maha <- knockdown:::.mahalanobis(x)
  ord.maha <- order(maha)
  testthat::expect_true(all(ord.x == ord.maha))
})

testthat::test_that("prioritization is computed correctly", {
  d <- data.table::data.table(A=letters[seq(10)], Readout=seq(10),
                              Pval=0, Qval=0)
  priorit <- knockdown:::.prioritize.chisq.statistic(d, 5, 1, 1)
  testthat::expect_true(all(priorit$Readout >= 5))
})

testthat::test_that("prioritization is computed correctly when Qval set", {
  d <- data.table::data.table(A=letters[seq(10)], Readout=seq(10),
                              Pval=0, Qval=seq(0, 1, length.out=10))
  priorit <- knockdown:::.prioritize.chisq.statistic(d, 5, 1, 0.5)
  testthat::expect_true(all(priorit$Qval <= .5))
})

testthat::test_that("chisq is correct", {
  c <- knockdown:::.chisq(1)
  testthat::expect_equal(0.3173105, c,  tolerance=.001)
})
