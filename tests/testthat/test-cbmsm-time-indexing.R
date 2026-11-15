test_that("CBMSM and balance.CBMSM work when time is not 1-based", {
  data(Blackwell)
  bw_shifted <- Blackwell
  bw_shifted$time <- bw_shifted$time + 1999

  form0 <- "d.gone.neg ~ d.gone.neg.l1 + camp.length"
  fit_base <- CBMSM(
    formula = form0, time = Blackwell$time, id = Blackwell$demName,
    data = Blackwell, type = "MSM", iterations = NULL, twostep = TRUE,
    msm.variance = "approx", time.vary = FALSE, init = "glm"
  )
  fit_shifted <- CBMSM(
    formula = form0, time = bw_shifted$time, id = bw_shifted$demName,
    data = bw_shifted, type = "MSM", iterations = NULL, twostep = TRUE,
    msm.variance = "approx", time.vary = FALSE, init = "glm"
  )

  expect_equal(unname(fit_shifted$weights), unname(fit_base$weights))
  expect_equal(unname(fit_shifted$fitted.values), unname(fit_base$fitted.values))
  expect_equal(unname(fit_shifted$glm.weights), unname(fit_base$glm.weights))

  bal_base <- balance.CBMSM(fit_base)
  bal_shifted <- balance.CBMSM(fit_shifted)
  expect_equal(bal_shifted$Balanced, bal_base$Balanced)
})
