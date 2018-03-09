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


testthat::test_that("hm object prints", {
  testthat::expect_output(show(hm.fit))
})

testthat::test_that("hm object plots", {
  suppressWarnings(s <- plot(hm.fit))
  testthat::expect_silent(s[[1]])
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
  pars <- params(fit)$formula
  print(pars)
  testthat::expect_true(is.character(pars))
  testthat::expect_true(pars == "Readout ~ Condition + (1 | GeneSymbol) + (1 | Condition:GeneSymbol) + (1 | ScreenType) + (1 | Condition:ScreenType)")
})


testthat::test_that("hm object returns params", {
  testthat::expect_true(class(params(hm.fit))[1] == "list")
})


testthat::test_that("hm object returns modelfit", {
  testthat::expect_true(class(modelFit(hm.fit))[1] == "list")
})


testthat::test_that("hm object returns data", {
  testthat::expect_true(class(dataSet(hm.fit))[1] == "data.table")
})


testthat::test_that("hm object returns gene effects", {
  testthat::expect_true(class(geneEffects(hm.fit))[1] == "data.table")
})


testthat::test_that("hm object returns gene hits", {
  testthat::expect_true(class(geneHits(hm.fit))[1] == "data.table")
})


testthat::test_that("hm object returns nested gene effects", {
  testthat::expect_true(class(nestedGeneEffects(hm.fit))[1] == "data.table")
})


testthat::test_that("hm object returns nested gene hits", {
  testthat::expect_true(class(nestedGeneHits(hm.fit))[1] == "data.table")
})


testthat::test_that("hm object returns bootstrapping boolean", {
  testthat::expect_false(isBootstrapped(hm.fit))
})


testthat::test_that("hm object returns inference", {
  testthat::expect_true(inference(hm.fit) ==
                        perturbatr:::inferenceTypes()$MIXED.MODEL)
})


testthat::test_that("hm bootstrapping throws with boostrap < 10", {
  testthat::expect_error(
    hm(rnaiscreen, effect.size=0.01, qval.threshold=1, bootstrap.cnt=4))
})


testthat::test_that("hm bootstrapping works", {
  ft <- hm(rnaiscreen, effect.size=0.01, qval.threshold=1, bootstrap.cnt=10)
  testthat::expect_true(isBootstrapped(ft))
  testthat::expect_true(all(!is.na(geneHits(ft))))
})
