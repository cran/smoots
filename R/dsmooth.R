#' Data-driven Local Polynomial for the Trend's Derivatives in Equidistant Time Series
#'
#'  This function runs through an iterative process in order to find the
#'  optimal bandwidth for the nonparametric estimation of the first or second
#'  derivative of the trend in an equidistant time series (with short-memory
#'  errors) and subsequently employs the obtained bandwidth via local
#'  polynomial regression.
#'
#'@param y a numeric vector that contains the time series ordered from past to present.
#'@param d an integer 1 or 2 that defines the order of derivative.
#'@param mu an integer 0, ..., 3 that represents the smoothness parameter of
#'the kernel weighting function and thus defines the kernel function that will
#'be used within the local polynomial regression; is set to \emph{1} by
#'default.
#'
#'\tabular{cl}{
#'\strong{Number} \tab \strong{Kernel}\cr
#'\emph{0} \tab Uniform Kernel\cr
#'\emph{1} \tab Epanechnikov Kernel\cr
#'\emph{2} \tab Bisquare Kernel\cr
#'\emph{3} \tab Triweight Kernel
#'}
#'@param pp an integer 1 (local linear regression) or 3 (local cubic
#'regression) that indicates the order of polynomial upon which c_[f], i.e. the
#'variance factor, will be calculated.
#'@param bStart.p a numeric object that indicates the starting value of the
#'  bandwidth for the iterative process for the calculation of c_[f]; should be
#'  0 < bStart.p < 0.5; is set to \emph{0.15} by default.
#'@param bStart a numeric object that indicates the starting value of the
#'  bandwidth for the iterative process; should be 0 < bStart < 0.5; is set to
#'  \emph{0.15} by default.
#'
#'@details
#'The iterative-plug-in (IPI) algorithm, which numerically minimizes the
#'Asymptotic Mean Squared Error (AMISE), was proposed by Feng et al. (2019).
#'
#'Define I[m^(k)] = int_[c_[b]]^[d_[b]] [m^(k)(x)]^2 dx,
#'
#'beta_[v,k] =  int_[-1]^[1] u^k K(u)du and
#'
#'R(K) = int_[-1]^[1] K^2(u)du.
#'
#'The AMISE is then
#'
#'AMISE(h) = h^(2(k-v)) * ( I[m^(k)]beta^2 / [k]^2 )
#'         + ( 2pi * c_[f](d_[b]-c_[b])R(K) / nh^(2v+1) )
#'
#'with h being the bandwidth, k = p + 1 being the order of the equivalent
#'kernel, v being the order of derivative, 0 <= c_[b] < d_[b] <= 1, n being the
#'number of observations, c_[f] being the variance factor and K_[(v,k)](u)
#'being the k-th order equivalent kernel obtained for the estimation of m^[(v)]
#'in the interior. m^[(v)] is the v-th order derivative (v = 0, 1, 2, ...) of
#'the nonparametric trend.
#'
#'The variance factor c_[f] is first obtained from a pilot-estimation of the
#'time series' nonparametric trend (v = 0) with polynomial order p_[p].
#'The estimate is then plugged into the iterative procedure for estimating
#'the first or second derivative (v = 1 or v = 2). For further details on the
#'asymptotic theory or the algorithm, we refer the user to Feng, Fritz and
#'Gries (2019) and Feng et al. (2019).
#'
#'The function itself is applicable in the following way: Based on a data input
#'\emph{y}, an order of polynomial \emph{pp} for the variance factor estimation
#'procedure, a starting value for the relative bandwidth \emph{bStart.p} in the
#'variance factor estimation procedure, a kernel function defined by the
#'smoothness parameter \emph{mu} and a starting value for the relative
#'bandwidth \emph{bStart} in the bandwidth estimation procedure, an optimal
#'bandwidth is numerically calculated for the derivative of order \emph{d}.
#'In fact, aside from the input vector \emph{y}, every argument has a default
#'setting that can be adjusted for the individual case. However, it is
#'recommended to initially use the default values for the estimation of the
#'first derivative and adjust the argument \emph{d} to \emph{d = 2} for the
#'estimation of the second derivative. Following Feng, Gries and Fritz (2019),
#'the initial bandwidth does not affect the resulting optimal bandwidth in
#'theory. However in practice, local minima of the AMISE can influence the
#'results. Therefore, the default starting bandwidth is set to 0.15, the
#'suggested starting bandwidth by Feng, Gries and Fritz (2019) for the
#'data-driven  estimation of the first derivative. The recommended initial
#'bandwidth for the second derivative, however, is 0.2 and not 0.15. Thus, if
#'the algorithm does not give suitable results (especially for \emph{d = 2}),
#'the adjustment of the initial bandwidth might be a good starting point.
#'Analogously, the default starting bandwith for the trend estimation for the
#'variance factor is \emph{bStart.p = 0.15}, although according to Feng, Gries
#'and Fritz (2019), \emph{bStart.p = 0.1} is suggested for \emph{pp = 1} and
#'\emph{bStart.p = 0.2} for \emph{pp = 3}. The default is therefore a
#'compromise between the two suggested values. For more specific information
#'on the input arguments consult the section \emph{Arguments}.
#'
#'After the bandwidth estimation, the nonparametric derivative of the series
#'is calulated with respect to the obtained optimal bandwidth by means of a
#'local polynomial regression. The output object is then a list that contains,
#'among other components, the original time series, the estimates of the
#'derivative and the estimated optimal bandwidth.
#'
#'The default print method for this function delivers key numbers such as
#'the iteration steps and the generated optimal bandwidth rounded to the fourth
#'decimal. The exact numbers and results such as the estimated nonparametric
#'trend series are saved within the output object and can be addressed via the
#'\emph{$} sign.
#'
#'@return The function returns a list with different components:
#'
#'\describe{
#'\item{b0}{the optimal bandwidth chosen by the IPI-algorithm.}
#'\item{bStart}{the starting bandwidth for the local polynomial
#'regression based derivative estimation procedure; input argument.}
#'\item{bStart.p}{the starting bandwidth for the nonparametric trend estimation
#'that leads to the variance factor estimate; input argument.}
#'\item{cf0}{the estimated variance factor.}
#'\item{InfR}{the inflation rate setting.}
#'\item{iterations}{the bandwidths of the single iterations steps}
#'\item{mu}{the smoothness parameter of the second order kernel; input
#'argument.}
#'\item{n}{the number of observations.}
#'\item{niterations}{the total number of iterations until convergence.}
#'\item{orig}{the original input series; input argument.}
#'\item{p}{the order of polynomial for the local polynomial
#'regression used within derivative estimation procedure.}
#'\item{pp}{the order of polynomial for the local polynomial
#'regression used in the variance factor estimation; input argument.}
#'\item{res}{the estimated residual series.}
#'\item{ws}{the weighting systems used within the local polynomial regression;
#'only exists, if the final smoothing is done via a local polynomial
#'regression.}
#'\item{ye}{the nonparametric estimates of the derivative.}
#'}
#'
#'@export
#'
#'@references
#' Feng, Y., Gries, T. and Fritz, M. (2019). Data-driven
#' local polynomial for the trend and its derivatives in economic time
#' series. Discussion Paper. Paderborn University.
#'
#' Feng, Y., Gries, T., Letmathe, S., and Schulz, D. (2019). The smoots package
#' in R for semiparametric modeling of trend stationary time series. Discussion
#' Paper. Paderborn University.
#'
#'@author
#'\itemize{
#'\item Yuanhua Feng (Department of Economics, Paderborn University), \cr
#'Author of the Algorithms \cr
#'Website: \url{https://wiwi.uni-paderborn.de/en/dep4/feng/}
#'\item Dominik Schulz (Student Assistant) (Department of Economics, Paderborn
#'University), \cr
#'Package Creator and Maintainer
#'}
#'
#'@examples
#'# Logarithm of test data
#'test_data <- gdpUS
#'y <- log(test_data$GDP)
#'t <- seq(from = 1947, to = 2019.25, by = 0.25)
#'
#'# Applied dsmooth function for the trend's first derivative
#'result_d <- dsmooth(y, d = 1, mu = 1, pp = 1, bStart.p = 0.1, bStart = 0.15)
#'estim <- result_d$ye
#'
#'# Plot of the results
#'plot(t, estim, xlab = "Year", ylab = "First derivative", type = "l",
#'  main = "Estimated first derivative of the trend for log-quarterly US-GDP, Q1 1947 - Q2 2019",
#'  cex.axis = 0.8, cex.main = 0.8, cex.lab = 0.8, bty = "n")
#'
#'# Print result
#'result_d


# The main function "dsmooth"------------------------------------------

# d = 1 or 2 indicate iterative plug-in for m' or m'', respectively
dsmooth <- function(y, d = c(1, 2), mu = c(0, 1, 2, 3), pp = c(1, 3),
                    bStart.p = 0.15, bStart = 0.15) {

  if (all(d == c(1, 2))) d <- 1
  if (all(mu == c(0, 1, 2, 3))) mu <- 1
  if (all(pp == c(1, 3))) pp <- 1

  if (!(d %in% c(1, 2))) {
    stop("Input of argument 'd' incorrect. It must be set to either 1 or 2.")
  }
  if (!(pp %in% c(1, 3))) {
    stop("Input of argument 'pp' incorrect. It must be set to either 1 or 3.")
  }
  if (!(is.numeric(bStart.p) && length(bStart.p) == 1)) {
    stop("Input of argument 'bStart.p' incorrect.
            One numerical value between 0 and 0.5 is needed.")
  }
  if (bStart.p <= 0 || bStart.p >= 0.5) {
    message("NOTE: 'bStart.p' between 0 and 0.5 is recommended.")
  }
  if (!(is.numeric(bStart) && length(bStart) == 1)) {
    stop("Input of argument 'bStart' incorrect.
         One numerical value between 0 and 0.5 is needed.")
  }
  if (bStart <= 0 || bStart >= 0.5) {
    message("NOTE: 'bStart' between 0 and 0.5 is recommended.")
  }

  if (pp == 1) {result.p <- msmooth(y, 1, mu, bStart.p, "A")}
  if (pp == 3) {result.p <- msmooth(y, 3, mu, bStart.p, "B")}
  cf0 <- result.p$cf0

  # Using the knn idea, bb=1, or not
  bb <- 1  # default
  if (d == 1) {InfR <- "Nai"} else {InfR <- "Var"}  # default
  cb <- 0.05
  # In this code ccf and bStart will be provided by the user
  # Input parameters
  n <- length(y)
  p <- d + 1
  k <- p + 1
  pd <- p + 2
  runc <- 1
  n1 <- trunc(n * cb)

  # New method for the kernel constants with p = 1 or 3------------------------

  # Kernel constants-----------------------------------------------------------
  m <- 1000000  # for the numerical integral
  u <- (-m:m) / (m + 0.5)
  # For d = 1, the four 3rd-order kernels for m' (Table 5.7, Mueller, 1988)

  wkp <- lookup$lookup_3kerns[(mu + 1), d][[1]](u)

  Rp <- sum(wkp ** 2) / m
  mukp <- sum((u ** k) * wkp) / m

  # Two constant in the bandwidth
  c1 <- factorial(k) ** 2 * (2 * d + 1) / (2 * (k - d))  # for d = 1 or 2
  c2 <- (1 - 2 * cb) * Rp / (mukp) ** 2

  steps <- rep(NA, 40)

  # The main iteration---------------------------------------------------------

  noi <- 40  # maximal number of iterations

  for (i in 1:noi) {
    if (runc == 1) {
      if (i > 1) {bold1 <- bold}
      if (i == 1) {bold <- bStart} else {bold <- bopt}

      # Look up EIM inflation rates from the lookup table
      bd <- lookup$InfR2_lookup[d, InfR][[1]](bold)

      if (bd >= 0.49) {bd <- 0.49}
      yed <- gsmooth(y, k, pd, mu, bd, bb)$ye
      I2 <- sum(yed[(n1 + 1):(n - n1)] ** 2) / (n - 2 * n1)
	    c3 <- cf0 / I2

      if (d == 1) {
        bopt = (c1 * c2 * c3) ** (1 / 7) * n ** (-1 / 7)
        if (bopt < n ** (-7 / 9)) {bopt <- n ** (-7 / 9)}
      }
      if (d == 2) {
        bopt <- (c1 * c2 * c3) ** (1 / 9) * n ** (-1 / 9)
        if (bopt < n ** (-9 / 11)) {bopt <- n ** (-9 / 11)}
      }
      if (bopt > 0.49) {bopt <- 0.49}
      steps[i] = bopt

      if (i > 2 && abs(bold - bopt) / bopt < 1 / n) {runc <- 0}
      if (i > 3 && abs(bold1 - bopt) / bopt < 1 / n){
        bopt <- (bold + bopt) / 2
        runc <- 0
      }
    }
  }

  # Smooth with the selected bandwidth-----------------------------------------

  if (d == 1 && bopt < n ** (-7 / 9)) {bopt <- n ** (-7 / 9)}
  if (d == 2 && bopt < n ** (-9 / 11)) {bopt <- n ** (-9 / 11)}
  if (bopt > 0.49) {bopt <- 0.49}
  est.opt = gsmooth(y, d, p, mu, bopt, bb)
  ye <- est.opt$ye
  ws <- est.opt$ws
  # Final results
  results <- list(ye = ye, orig = y, b0 = bopt, ws = ws,
                  cf0 = cf0, deriv = d, n = length(y),
                  iterations = steps[!is.na(steps)],
                  niterations = length(steps[!is.na(steps)]),
                  pp = pp, bStart.p = bStart.p, bStart = bStart,
                  InfR = InfR, mu = mu, p = d + 1)

  class(results) <- "smoots"
  attr(results, "function") <- "dsmooth"

  drop(results)
}
# End of the function
