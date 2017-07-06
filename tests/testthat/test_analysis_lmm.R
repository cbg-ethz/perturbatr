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


context("analysis")

data(rnaiscreen)
rnai.norm <- preprocess(rnaiscreen, normalize="z.score")

testthat::test_that("lmm model weights get set correctly", {
  mat <- set.lmm.model.data(rnai.norm, weights=list("pooled"=2, "single"=1))
  testthat::expect_true(
    all(mat@.data$Weight[mat@.data$Design == "pooled"] == 2))
})

testthat::test_that("lmm model weights throws", {
  testthat::expect_error(
    set.lmm.model.data(rnai.norm, weights=list("pooled"=2, "unpooled"=1)))
})

testthat::test_that("lmm model group is set correctly", {
  mat <- set.lmm.model.data(rnai.norm, weights=list("pooled"=2, "single"=1))
  testthat::expect_equal(as.character(mat@.data$VG),
                         paste(sep=":", mat@.data$Virus, mat@.data$GeneSymbol))
})
