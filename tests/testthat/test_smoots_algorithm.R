Xt <- log(smoots::gdpUS$GDP)
est1 <- smoots::msmooth(Xt)
est2 <- smoots::dsmooth(Xt)
est3 <- smoots::dsmooth(Xt, d = 2)
est4 <- smoots::tsmooth(Xt, p = 3, InfR = "Nai", mu = 3, bb = 0, Mcf = "AR",
  bvc = "N")
est5 <- smoots::dsmooth(Xt, d = 2, mu = 3, pp = 3)
b1 <- est1$b0
b2 <- est2$b0
b3 <- est3$b0
b4 <- est4$b0
b5 <- est5$b0
ye1.head <- head(est1$ye, 3)
ye1.mid <- est1$ye[151:153]
ye1.tail <- tail(est1$ye, 3)
ye2.head <- head(est2$ye, 3)
ye2.mid <- est2$ye[151:153]
ye2.tail <- tail(est2$ye, 3)
ye3.head <- head(est3$ye, 3)
ye3.mid <- est3$ye[151:153]
ye3.tail <- tail(est3$ye, 3)
ye4.head <- head(est4$ye, 3)
ye4.mid <- est4$ye[151:153]
ye4.tail <- tail(est4$ye, 3)
ye5.head <- head(est5$ye, 3)
ye5.mid <- est5$ye[151:153]
ye5.tail <- tail(est5$ye, 3)

b1.test <- 0.132322665173448
b2.test <- 0.205425095576514
b3.test <- 0.254715050452228
b4.test <- 0.100171310780351
b5.test <- 0.267922433080052
ye1.head.test <- c(7.61876243092676,
                   7.62807830529724,
                   7.63739258039812)
ye1.mid.test <- c(8.94726967204048,
                  8.95501578931484,
                  8.96272389387410)
ye1.tail.test <- c(9.82319520210691,
                   9.82770154886513,
                   9.83220887427568)
ye2.head.test <- c(2.73795232321980,
                   2.73875814816755,
                   2.73956791840870)
ye2.mid.test <- c(2.22043267682466,
                  2.22321167556111,
                  2.22612401451022)
ye2.tail.test <- c(0.777043594443846,
                   0.760964141212400,
                   0.744898872849721)
ye3.head.test <- c(4.87055549064138,
                   4.79230656799487,
                   4.71391958703809)
ye3.mid.test <- c(0.863413946823043,
                  0.875282050169084,
                  0.884475520161740)
ye3.tail.test <- c(-4.06446240594676,
                   -4.04268398885208,
                   -4.02072357326614)
ye4.head.test <- c(7.62203280122757,
                   7.62225187487096,
                   7.62257278971999)
ye4.mid.test <- c(8.93100433043105,
                  8.94197739706743,
                  8.95320516281173)
ye4.tail.test <- c(9.84123927973195,
                   9.84859560776920,
                   9.85628471403288)
ye5.head.test <- c(0.697781105931834,
                   0.737602242053834,
                   0.776602464045702)
ye5.mid.test <- c(0.755528599113370,
                  0.755560697533622,
                  0.755479249607318)
ye5.tail.test <- c(2.00662550738630,
                   2.18262038279371,
                   2.35943401286232)

test_that("Trend estimation consistency via msmooth", {
  expect_equal(b1, b1.test)
  expect_equal(ye1.head, ye1.head.test)
  expect_equal(ye1.mid, ye1.mid.test)
  expect_equal(ye1.tail, ye1.tail.test)
})

test_that("First derivative estimation consistency via dsmooth", {
  expect_equal(b2, b2.test)
  expect_equal(ye2.head, ye2.head.test)
  expect_equal(ye2.mid, ye2.mid.test)
  expect_equal(ye2.tail, ye2.tail.test)
})

test_that("Second derivative estimation consistency via dsmooth", {
  expect_equal(b3, b3.test)
  expect_equal(ye3.head, ye3.head.test)
  expect_equal(ye3.mid, ye3.mid.test)
  expect_equal(ye3.tail, ye3.tail.test)
})

test_that("Trend estimation consistency via tsmooth with more complex settings", {
  expect_equal(b4, b4.test)
  expect_equal(ye4.head, ye4.head.test)
  expect_equal(ye4.mid, ye4.mid.test)
  expect_equal(ye4.tail, ye4.tail.test)
})

test_that("Second derivative estimation consistency via dsmooth with more complex settings", {
  expect_equal(b5, b5.test)
  expect_equal(ye5.head, ye5.head.test)
  expect_equal(ye5.mid, ye5.mid.test)
  expect_equal(ye5.tail, ye5.tail.test)
})
