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


test.file <- system.file("extdata", "testscreen.rds", package = "perturbatr")
rnaiscreen <- readRDS(test.file)
hm.fit  <- hm(rnaiscreen)


testthat::test_that("hm object prints", {
  testthat::expect_output(show(hm.fit))
})


testthat::test_that("hm object plots", {
  suppressWarnings(s <- plot(hm.fit))
  testthat::expect_silent(s[[1]])
})


testthat::test_that("hm has correct inference", {
  fit  <- hm(rnaiscreen,
             formula=Readout ~ Condition+(1|GeneSymbol)+(1|Condition:GeneSymbol) + (1|ScreenType)+(1|Condition:ScreenType))
  pars <- params(fit)$formula
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
  testthat::expect_true(class(dataSet(hm.fit))[1] == "tbl_df")
})


testthat::test_that("hm object returns gene effects", {
  testthat::expect_true(class(geneEffects(hm.fit))[1] == "tbl_df")
})


testthat::test_that("hm object returns nested gene effects", {
  testthat::expect_true(class(nestedGeneEffects(hm.fit))[1] == "tbl_df")
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
    hm(rnaiscreen, bootstrap.cnt=4))
})


testthat::test_that("hm bootstrapping works", {
  ft <- hm(rnaiscreen, bootstrap.cnt=10)
  testthat::expect_true(isBootstrapped(ft))
  testthat::expect_true(all(!is.na(geneEffects(ft))))
})


testthat::test_that("hm throws with no gene:condition random effect", {
  testthat::expect_error(
    hm(rnaiscreen,
       formula = Readout ~ Condition+(1|GeneSymbol)+ (1|ScreenType)+(1|Condition:ScreenType)))
})


testthat::test_that("hm throws with no gene random effect", {
  testthat::expect_error(
    hm(rnaiscreen,
       formula = Readout ~ Condition+(1|Condition:GeneSymbol) + (1|ScreenType)+(1|Condition:ScreenType)))
})


testthat::test_that("hm throws with no condition fixed effect", {
  testthat::expect_error(
    hm(rnaiscreen,
       formula = Readout ~ (1|GeneSymbol)+(1|Condition:GeneSymbol) + (1|ScreenType)+(1|Condition:ScreenType)))
})


testthat::test_that("hm throws with completely false model", {
  testthat::expect_error(
    hm(rnaiscreen,
       formula = (1|ScreenType)+(1|Condition:ScreenType)))
})
