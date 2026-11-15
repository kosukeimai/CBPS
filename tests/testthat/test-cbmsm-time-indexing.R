test_that("CBMSM and balance.CBMSM do not depend on how time is coded", {
  data(Blackwell)
  form0 <- "d.gone.neg ~ d.gone.neg.l1 + camp.length"
  fit_with_time <- function(time) {
    d <- Blackwell
    d$time <- time
    set.seed(12345)
    CBMSM(
      formula = form0, time = d$time, id = d$demName,
      data = d, type = "MSM", iterations = NULL, twostep = TRUE,
      msm.variance = "approx", time.vary = FALSE, init = "glm"
    )
  }

  t0 <- Blackwell$time
  fit_base <- fit_with_time(t0)
  bal_base <- balance(fit_base)

  codings <- list(
    "calendar years" = t0 + 1999,
    "0-based" = t0 - 1,
    "gaps" = 2 * t0 - 1
  )
  for (nm in names(codings)) {
    fit <- fit_with_time(codings[[nm]])
    expect_equal(unname(fit$weights), unname(fit_base$weights), info = nm)
    expect_equal(unname(fit$fitted.values), unname(fit_base$fitted.values), info = nm)
    expect_equal(unname(fit$glm.weights), unname(fit_base$glm.weights), info = nm)
    expect_equal(unname(fit$treat.hist), unname(fit_base$treat.hist), info = nm)
    expect_equal(sort(unique(fit$time)), sort(unique(codings[[nm]])), info = nm)
    expect_equal(balance(fit)$Balanced, bal_base$Balanced, info = nm)
  }
})
