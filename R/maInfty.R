maInfty <- function(beta, alpha, m = 1000) {
  p <- length(beta)
  q <- length(alpha)
  if (m - q > 0) {
    times <- m - q
  } else {
    times <- 0
  }
  alpha.s <- c(alpha, rep(0, times = times))
  c.out <- c(rep(0, times = p - 1), 1, rep(NA, times = m))
  lc <- length(c.out)

  if (m >= 1) {
    for (i in (p + 1):(m + p)) {
      c.out[i] <- beta %*% c.out[(i - 1):(i - p)] + alpha.s[i - p]
    }
  }

  c.out[p:lc]
}
