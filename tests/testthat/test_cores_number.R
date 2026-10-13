test_that("Ensure that number of cores is handled correctly in bootstrap", {
  est <- msmooth(log(gdpUS$GDP))
  # Error message for cores < 1
  expect_error({
    fc <- modelCast(est, p = 1, method = "boot", cores = 0, pb = FALSE)
  }, "The argument 'cores' must be a single integer value > 0 or NULL.")
  # No errors for values >= 1
  expect_no_error({
    fc <- modelCast(est, p = 1, method = "boot", cores = 1, pb = FALSE)
  })
  expect_no_error({
    fc <- modelCast(est, p = 1, method = "boot", cores = 2, pb = FALSE)
  })
  # No errors on default
  expect_no_error({
    fc <- modelCast(est, p = 1, method = "boot", pb = FALSE)
  })
})
