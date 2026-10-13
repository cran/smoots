#'
#'@importFrom stats arima
#'

# Function for selecting M and estimating cf nonparametrically-----------------

# removed and now written in C++

#------------------------------------------------------------------------------


# Function for estimating cf0 from AR(0) to AR(p.max), e.g. p.max = 5----------

cf0.AR.est <- function(Xt, nar_min, nar_max, arma_mean) {
  p.min <- nar_min
  p.max <- nar_max
  p.sel <- p.min:p.max
  BIC <- (p.sel) * 0
  n = length(Xt)
  for (l in p.sel) {
    AR.l <- suppressWarnings(arima(Xt, order = c(l, 0, 0), include.mean = arma_mean))
    BIC[l - p.min + 1] <- -2 * AR.l$loglik + log(n) * (l + 2) - log(n) * 2
  }
  p.BIC <- p.sel[BIC == min(BIC)]
  AR.BIC <- suppressWarnings(arima(Xt, order = c(p.BIC, 0, 0), include.mean = arma_mean))
  sc.AR <- 0
  if (p.BIC > 0) {sc.AR <- sum(AR.BIC$coef[1:p.BIC])}
  cf0.AR <- (1 / (1 - sc.AR)) ** 2 * AR.BIC$sigma2

  results <- list(cf0.AR = cf0.AR, AR.BIC = AR.BIC, p.BIC = p.BIC)

  results
}
#### end of the cf0.AR.est function

### Estimation of optimal ARMA model and the resulting variance factor

cf0.ARMA.est <- function(Xt, nar_min, nar_max, nma_min, nma_max, arma_mean) {
  p.min <- nar_min
  p.max <- nar_max
  q.min <- nma_min
  q.max <- nma_max
  p.arma <- rep(p.min:p.max, times = length(q.min:q.max))
  q.arma <- rep(q.min:q.max, each = length(p.min:p.max))
  bic <- mapply(function(x, p0, q0) {
    arma.nobs <- length(x)
    arma.est <- suppressWarnings(arima(x, order = c(p0, 0, q0), include.mean = arma_mean))
    -2 * arma.est$loglik + log(arma.nobs) * (p0 + q0)
  }, p0 = p.arma, q0 = q.arma,
  MoreArgs = list(x = Xt))
  opt <- which(bic == min(bic))
  p.BIC <- p.arma[opt]
  q.BIC <- q.arma[opt]
  ARMA.BIC <- suppressWarnings(arima(Xt, order = c(p.BIC, 0, q.BIC), include.mean = arma_mean))
  if (p.BIC != 0) {
    sc.AR <- sum(ARMA.BIC$coef[1:p.BIC])
  } else {
    sc.AR <- 0
  }
  if (q.BIC != 0) {
    sc.MA <- sum(ARMA.BIC$coef[(p.BIC + 1):(p.BIC + q.BIC)])
  } else {
    sc.MA <- 0
  }
  cf0.ARMA <- ((1 + sc.MA) / (1 - sc.AR)) ^ 2 * ARMA.BIC$sigma2
  results <- list(cf0.ARMA = cf0.ARMA, ARMA.BIC <- ARMA.BIC,
                  p.BIC = p.BIC, q.BIC = q.BIC)

  results
}

### Estimation of the optimal MA model and the resulting variance factor

cf0.MA.est <- function(Xt, nma_min, nma_max, arma_mean) {
  q.min <- nma_min
  q.max <- nma_max
  q.arma <- q.min:q.max
  bic <- sapply(q.arma, FUN = function(x, q0) {
    ma.nobs <- length(x)
    ma.est <- suppressWarnings(arima(x, order = c(0, 0, q0), include.mean = arma_mean))
    -2 * ma.est$loglik + log(ma.nobs) * q0
  }, x = Xt)
  opt <- which(bic == min(bic))
  q.BIC <- q.arma[opt]
  MA.BIC <- suppressWarnings(arima(Xt, order = c(0, 0, q.BIC), include.mean = arma_mean))
  if (q.BIC != 0) {
    sc.MA <- sum(MA.BIC$coef[1:q.BIC])
  } else {
    sc.MA <- 0
  }
  cf0.MA <- (1 + sc.MA) ^ 2 * MA.BIC$sigma2
  results <- list(cf0.MA = cf0.MA, MA.BIC <- MA.BIC,
                  q.BIC = q.BIC)

  results
}
