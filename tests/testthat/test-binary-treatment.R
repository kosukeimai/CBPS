test_that("binary fits retain the observation names and requested treatment", {
  set.seed(104)
  x <- rnorm(120)
  z <- rnorm(120)
  active <- rbinom(120, 1, plogis(x / 3 - z / 4))
  group <- factor(ifelse(active == 1, "A", "B"), levels = c("B", "A"))
  ids <- paste0("unit-", seq_along(x))
  names(group) <- ids
  design <- cbind(1, x, z)

  for (mode in c("ate", "att2", "optimal")) {
    optimal <- mode == "optimal"
    reference <- if (mode == "att2") 1 - active else active
    fit <- CBPS.fit(
      treat = group, X = design,
      baselineX = if (optimal) as.matrix(x), diffX = if (optimal) as.matrix(z),
      ATT = if (mode == "att2") 2 else 0, method = "over",
      iterations = 1000, standardize = TRUE, twostep = TRUE,
      sample.weights = rep(1, length(x))
    )
    expect_identical(fit$y, setNames(as.numeric(reference), ids))
    p <- plogis(drop(design %*% coef(fit)))
    raw <- if (mode == "att2") {
      ifelse(reference == 1, 1, p / (1 - p))
    } else {
      ifelse(reference == 1, 1 / p, 1 / (1 - p))
    }
    w <- raw / ave(raw, reference, FUN = sum)
    expect_equal(unname(fit$fitted.values), unname(p), tolerance = 1e-12)
    expect_equal(unname(fit$weights), unname(w), tolerance = 1e-12)
  }
})

test_that("unnamed binary character codes return an unnamed treatment", {
  set.seed(105)
  design <- cbind(1, x = rnorm(100), z = rnorm(100))
  group <- rbinom(100, 1, 0.5)
  for (coded in list(as.character(group), ifelse(group == 1, "TRUE", "FALSE"))) {
    fit <- CBPS.fit(coded, design, NULL, NULL, 0, "over", 1000, TRUE, TRUE, rep(1, 100))
    expect_identical(fit$y, as.numeric(group))
  }
})
