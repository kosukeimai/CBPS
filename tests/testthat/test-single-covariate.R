test_that("CBPS preserves covariate names with a single covariate", {
  data(LaLonde)
  fit <- CBPS(treat ~ age, data = LaLonde, ATT = TRUE)

  expect_equal(rownames(fit$coefficients), c("(Intercept)", "age"))
  expect_equal(rownames(vcov(fit)), c("(Intercept)", "age"))
  expect_equal(colnames(vcov(fit)), c("(Intercept)", "age"))
  expect_equal(rownames(summary(fit)$coefficients), c("(Intercept)", "age"))
})
