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
# along with perturbatr If not, see <http://www.gnu.org/licenses/>.


context("clazz")

data(rnaiscreen)
test.dat <- dataSet(rnaiscreen) %>%
  dplyr::select(Condition, Replicate, Plate, RowIdx, ColIdx,
                GeneSymbol, Readout, Control, Perturbation) %>%
  as.data.frame

dat <- methods::as(test.dat, "PerturbationData")
dat.lmm <- hm(dat, effect.size=0.01)


testthat::test_that("data object has correct class", {
  dat  <- data.table(A=rnorm(10))
  dat2 <- data.frame(Condition=rnorm(10), Replicate=1, GeneSymbol="A")
  testthat::expect_error(methods::as(dat, "PerturbationData"))
  testthat::expect_error(methods::as(dat2, "PerturbationData"))
})


testthat::test_that("raw object has correct class", {
  testthat::expect_s4_class(dat, "PerturbationData")
})


testthat::test_that("raw object prints", {
  testthat::expect_output(show(rnaiscreen))
})


testthat::test_that("hm analysed object has correct class", {
  testthat::expect_s4_class(dat.lmm, "HMAnalysedPerturbationData")
})


testthat::test_that("hm analysed object has correct class", {
  testthat::expect_error(methods::new("AbstractAnalysedPerturbationData"))
})
