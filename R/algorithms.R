tsmooth.Alg <- function(y, p, mu, bStart, alg, method, nar_lim, nma_lim, arma_mean) {
  results <- tsmoothCalc(y, p, mu, lookup$lookup_alg["Mcf", alg],
                     lookup$lookup_alg["InfR", alg],
                     bStart, lookup$lookup_alg["bvc", alg], 1, 0.05, method,
                     nar_lim = nar_lim, nma_lim = nma_lim, arma_mean = arma_mean)
  results
}
