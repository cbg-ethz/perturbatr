# perturbationmodels: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2018 Simon Dirmeier
#
# This file is part of perturbationmodels
#
# perturbationmodels is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# perturbationmodels is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with perturbationmodels. If not, see <http://www.gnu.org/licenses/>.


context("hm analysis")


data(rnaiscreen)
rnai.norm <- preprocess(rnaiscreen, normalize="z.score")
hm.fit    <- hm(rnai.norm, effect.size=0.01, qval.threshold=1)

testthat::test_that("hm model weights get set correctly", {
  mat <- set.hm.model.data(rnai.norm, weights=list("pooled"=2, "single"=1))
  testthat::expect_true(
    all(mat@.data$Weight[mat@.data$Design == "pooled"] == 2))
})

testthat::test_that("hm model weights throws", {
  testthat::expect_error(
    set.hm.model.data(rnai.norm, weights=list("pooled"=2, "unpooled"=1)))
})

testthat::test_that("hm model group is set correctly", {
  mat <- set.hm.model.data(rnai.norm, weights=list("pooled"=2, "single"=1))
  testthat::expect_equal(as.character(mat@.data$VG),
                         paste(sep=":", mat@.data$Condition, mat@.data$GeneSymbol))
})

testthat::test_that("hm has correct inference", {
  testthat::expect_equal(hm.fit@.inference,
                         perturbatr:::.inference.types()$MIXED.MODEL)
})

testthat::test_that("hm has not been bootstrapped", {
  testthat::expect_false(hm.fit@.is.bootstrapped)
})

testthat::test_that("hm uses qval threshold of 1", {
  testthat::expect_equal(hm.fit@.params$qval.threshold, 1)
})

testthat::test_that("hm uses effect size threshold of .01", {
  testthat::expect_equal(hm.fit@.params$effect.size, .01, tolerance=0.001)
})
