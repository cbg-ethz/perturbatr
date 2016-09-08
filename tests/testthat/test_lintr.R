context("lints")

test_that("package has style", {
  lintr::expect_lint_free()
})
