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

testthat::test_that("warning is thrown when only two replicates", {
  testthat::expect_warning(invisible(knockdown::.t.test("gene", c(1, 2), 0)))
})

testthat::test_that("p value is one for two replicates", {
  options(warn=0)
  testthat::expect_warning(
    testthat::expect_equal(knockdown:::.t.test("gene", c(1, 2), 0)$p.value, 1)
  )
})

testthat::test_that("hit ratio for prioritization for t test is working", {
  dat <- rnaiscreen@.data %>% .[1:5] %>%
    dplyr::mutate(Pval = .01, Qval= .01, Readout = 10, GeneSymbol="a", Entrez=1)
  pri <- knockdown:::.prioritize.tstatistic(dat)
  testthat::expect_true(all(unlist(pri$HitRatio) == 1))
})

testthat::test_that("pval for prioritization for t test is working", {
  dat <- rnaiscreen@.data %>% .[1:5] %>%
    dplyr::mutate(Pval = .01, Qval= .01, Readout = 10, GeneSymbol="a", Entrez=1)
  pri <- knockdown:::.prioritize.tstatistic(dat)
  testthat::expect_equal(unlist(pri$Pval), metap::sumlog(rep(.01, 5))$p, tolerance=.001)
})
