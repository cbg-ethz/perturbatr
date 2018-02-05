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


context("analysis")


data(rnaiscreen)
rnai.norm <- preprocess(rnaiscreen, normalize="z.score")

testthat::test_that("filtering works correctly", {
  frn <- pertubation::filter(rnaiscreen, Virus=="V1")
  testthat::expect_true(all(frn@.data$Virus == "V1"))
})

testthat::test_that("binding works correctly", {
  rnait <- pertubation::rbind(rnaiscreen, rnaiscreen)
  testthat::expect_equal(nrow(rnait@.data), 2 * nrow(rnaiscreen@.data))
})

testthat::test_that("selecting works correctly", {
  rnait <- pertubation::select(rnaiscreen, Virus)
  testthat::expect_equal(ncol(rnait), 1)
})
