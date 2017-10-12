rm(list=ls())
library(CBPS)
library(testthat)
context("tests CBPS")

accuracy <- 0.0001

# set random seed
set.seed(12345)

test_that("tests CBMS on the Lalonde data", {
  #Load the LaLonde data
  data(LaLonde)
  #' ## Estimate CBPS via logistic regression
  fit <- CBPS(treat ~ age + educ + re75 + re74 + I(re75==0) + I(re74==0), data = LaLonde, ATT = TRUE)
  x <- vcov(fit)

  expect_that(dim(x), is_equivalent_to(c(7,7)))
  expect_equal(x["age", "age"], 0.01208418, tolerance = accuracy)
  expect_equal(x["re74", "re75"], -0.03098431, tolerance = accuracy)
})  
