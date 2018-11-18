library(testthat)
context("Test R and Cpp code calculate the same BinaryPoissonModel log-likelihood")

library(ape)
library(BinaryPoissonUsingSPLITT)


N <- 100

tree <- rtree(N)
x <- sample(c(0, 1), size = N, replace = TRUE)


test_that(
  "BinaryPoissonModelLogLik == BinaryPoissonModelLogLikCpp",
  expect_equal(BinaryPoissonModelLogLik(x, tree, 1, 0.02, 0.8),
               BinaryPoissonModelLogLikCpp(x, tree, 1, 0.02, 0.8))
)
