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


context("diffusion")


data(rnaiscreen)
hm.fit  <- hm(rnaiscreen, effect.size=0.01, qval.threshold=1)

testthat::test_that("diffusion runs smoothly", {
  graph.file <- system.file("extdata", "graph_file.tsv",
                            package = "perturbatr")
  testthat::expect_silent(diffuse(hm.fit, path=graph.file, r=0.1))
})

testthat::test_that("diffusion runs smoothly", {
  graph.file <- system.file("extdata", "graph_file.tsv",
                            package = "perturbatr")
  res <- diffuse(hm.fit, path=graph.file, r=0.1)
  testthat::expect_true(class(graph(res))[1] == "igraph")
  testthat::expect_true(class(graph(res))[1] == "igraph")
})
