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


context("hm analysis")


data(rnaiscreen)
hm.fit  <- hm(rnaiscreen, effect.size=0.01, qval.threshold=1)


testthat::test_that("hm has correct inference", {
  testthat::expect_equal(inference(hm.fit),
                         perturbatr:::inferenceTypes()$MIXED.MODEL)
})


testthat::test_that("hm has correct inference", {
  testthat::expect_equal(inference(hm.fit),
                         perturbatr:::inferenceTypes()$MIXED.MODEL)
})


testthat::test_that("hm has not been bootstrapped", {
  testthat::expect_false(isBootstrapped(hm.fit))
})


testthat::test_that("hm uses qval threshold of 1", {
  testthat::expect_equal(params(hm.fit)$qval.threshold, 1)
})


testthat::test_that("hm uses effect size threshold of .01", {
  testthat::expect_equal(params(hm.fit)$effect.size, .01, tolerance=0.001)
})

testthat::test_that("hm has correct inference", {
  fit  <- hm(rnaiscreen, effect.size=0.01, qval.threshold=1,
             formula=Readout ~ Condition+(1|GeneSymbol)+(1|Condition:GeneSymbol) +
                                (1|ScreenType)+(1|Condition:ScreenType))
  testthat::expect_silent(1)
})
