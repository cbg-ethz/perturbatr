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
# along with perturbatr. If not, see <http://www.gnu.org/licenses/>.


context("class")


graph.file <- system.file("extdata", "graph_file.tsv", package="perturbatr")
data(rnaiscreen)

test.dat <- rnaiscreen@.data %>%
  dplyr::filter(Condition == "V1") %>%
  dplyr::select(Condition, Replicate, Plate, RowIdx, ColIdx,
                GeneSymbol, Readout, Control, Perturbation) %>%
  as.data.frame

dat        <- methods::as(test.dat, "perturbation.data")

dat.norm  <- preprocess(rnaiscreen, normalize="z.score")
dat.tana  <- tstatistic(dat.norm)
dat.hana  <- hyper.statistic(dat.norm)
dat.lmm   <- hm(dat.norm, effect.size=0.01)
dat.diff  <- diffuse(dat.lmm, node.start.count=1, path=graph.file)

testthat::test_that("raw object has correct class", {
  testthat::expect_s4_class(dat, "perturbation.raw.data")
})


testthat::test_that("normalized object has correct class", {
  testthat::expect_s4_class(dat.norm, "perturbation.normalized.data")
})


testthat::test_that("tstatistic analysed object has correct class", {
  testthat::expect_s4_class(dat.tana, "perturbation.tstatistic.analysed")
})


testthat::test_that("hyperstatistic analysed object has correct class", {
  testthat::expect_s4_class(dat.hana, "perturbation.hyper.analysed")
})


testthat::test_that("hm analysed object has correct class", {
  testthat::expect_s4_class(dat.lmm, "perturbation.hm.analysed")
})
