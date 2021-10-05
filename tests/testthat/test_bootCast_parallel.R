set.seed(123)
Xt1 <- stats::arima.sim(model = list(ar = 0.73, ma = -0.61), n = 200) + 3.1
bootfc1 <- smoots::bootCast(Xt1, p = 1, q = 1, h = 5, cores = 2, pb = FALSE,
  include.mean = TRUE, it = 1000, export.error = TRUE, plot = FALSE)

fc1 <- unname(c(bootfc1$fcast[1, ]))
err1 <- c(head(bootfc1$err, 3))

fc.test1 <- c(3.03303934096287, 3.05110662932237, 3.0626046443142, 3.06992197639646, 3.07457872356892)

err.test1 <- c(-0.184844853423316, 0.990023327646193, 0.378347673164732, -1.07861872693332, 1.02869459617643, -0.928473039415859, -0.265946302474376, -0.0577952466215566, -0.386239385251191, -0.277829399585893, -0.510163869792313, 1.96714008005794, -0.304055774248956, 0.84327292768902, -0.251304619616104)

test_that("bootCast is consistent (mean != 0, any number of cores)", {
  expect_equal(fc1, fc.test1, tolerance = 1e-03)
  expect_equal(err1, err.test1, tolerance = 1e-03)
})

