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


context("analysis")


data(rnaiscreen)
ft <- hm(rnaiscreen, effect.size=0.01)


testthat::test_that("dataSet for PerturbationData has correct columns", {
  tr <- sort(perturbatr:::.requiredDataCols()) %in%
        sort(colnames(dataSet(rnaiscreen)))
  testthat::expect_true(all(tr))
})


testthat::test_that("filtering works correctly", {
  testthat::expect_true(class(dataSet(ft))[1] == "data.table")
})


testthat::test_that("filtering works correctly", {
  frn <- perturbatr::filter(rnaiscreen, Condition == "V1")
  testthat::expect_true(all(dataSet(frn)$Condition == "V1"))
})


testthat::test_that("binding works correctly", {
  rnait <- rbind(rnaiscreen, rnaiscreen)
  testthat::expect_equal(nrow(dataSet(rnait)), 2 * nrow(dataSet(rnaiscreen)))
})
