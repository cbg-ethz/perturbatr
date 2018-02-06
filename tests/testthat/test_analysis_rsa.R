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

context("hypergeometric hypothesis tests")


data(rnaiscreen)
rnai.norm <- preprocess(rnaiscreen, normalize="z.score")

testthat::test_that("error is thrown when summarization is not logical", {
  testthat::expect_error(hyper.statistic(rnai.norm, do.summarization="no"))
})

testthat::test_that("error is thrown when summarization and leven is sirna", {
  testthat::expect_error(hyper.statistic(rnai.norm, do.summarization=TRUE,
                                         level="sirna"))
})

testthat::test_that("different levels give different results", {
  all.genes    <- rep(paste0("g", 1:5), 4)
  all.sirnas   <- rep(letters[1:10], 2)
  plates       <- rep(1, 20)
  rows <- cols <- 1:20
  readouts     <- rnorm(20)
  h.sirna <- perturbatr:::.hypertest(all.genes, all.sirnas,
                                   plates, rows, cols, readouts, "sirna")
  h.gene  <- perturbatr:::.hypertest(all.genes, all.sirnas,
                                   plates, rows, cols, readouts, "gene")
  testthat::expect_true(any(h.sirna != h.gene))
})

testthat::test_that("hypertest on group works", {
  ranks    <- c(1, 2, 3)
  all.genes <- letters
  v <- data.table(A=perturbatr:::.hypertest.grp(ranks, all.genes)) %>%
      tidyr::separate(A, c("Pval", "Hit"), sep="_")
  testthat::expect_equal(length(unique(v$Pval)), 1)
  testthat::expect_equal(length(unique(v$Hit)), 1)
})
