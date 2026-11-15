test_that("balance.CBMSM standardizes Unweighted means using baseline, not bal", {
  data(Blackwell)
  form0 <- "d.gone.neg ~ d.gone.neg.l1 + camp.length"
  fit <- CBMSM(
    formula = form0, time = Blackwell$time, id = Blackwell$demName,
    data = Blackwell, type = "MSM", iterations = NULL, twostep = TRUE,
    msm.variance = "approx", time.vary = FALSE, init = "glm"
  )

  b <- balance.CBMSM(fit)
  n_groups <- ncol(b$Unweighted) / 2
  t1 <- which(fit$time == min(fit$time))
  denom <- sd(fit$glm.weights * fit$x[t1, 2])

  expected_std <- b$Unweighted[1, seq_len(n_groups)] / denom
  actual_std <- b$Unweighted[1, n_groups + seq_len(n_groups)]
  expect_equal(unname(actual_std), unname(expected_std))
})
