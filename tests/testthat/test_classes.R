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


context("class")


graph.file <- system.file("extdata", "graph_file.tsv", package="perturbation")
data(rnaiscreen)
test.dat <- rnaiscreen@.data %>%
  dplyr::filter(Virus=="V1") %>%
  dplyr::select(Replicate, Plate, RowIdx, ColIdx,
                GeneSymbol, Readout, Control, siRNAIDs) %>%
  as.data.frame
dat        <- as(test.dat, "perturbation.data")
replicates <- replicates(dat)
plates     <- plates(dat)
quality    <- quality(dat)

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


testthat::test_that("lmm analysed object has correct class", {
  testthat::expect_s4_class(dat.lmm, "perturbation.lmm.analysed")
})


testthat::test_that("plates has correct class", {
  testthat::expect_s4_class(plates, "perturbation.plates")
})


testthat::test_that("replicates has correct class", {
  testthat::expect_s4_class(replicates, "perturbation.replicates")
})


testthat::test_that("quality has correct class", {
  testthat::expect_s4_class(quality, "perturbation.quality")
})
