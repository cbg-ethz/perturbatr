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


testthat::test_that("read graph returns igraph object", {
  fl <- system.file("extdata", "graph_file.tsv", package="perturbatr")
  gr <- perturbatr:::read.graph(path=fl, graph=NULL)
  testthat::expect_true(class(gr)[1] == "igraph")
})


testthat::test_that("check columns does its job", {
  testthat::expect_silent(check.columns(rnaiscreen@dataSet, "Condition"))
})


testthat::test_that("check columns throws on wrong cols", {
  testthat::expect_error(check.columns(rnaiscreen@dataSet, "wrong col"))
})


testthat::test_that("effect matrixes returns correct", {
  ef <- effect.matrices(ft)
  testthat::expect_true(class(ef) == "list")
  testthat::expect_true(class(ef$gene.effects)[1] == "data.table")
  testthat::expect_true(class(ef$nested.gene.effects)[1] == "data.table")
})
