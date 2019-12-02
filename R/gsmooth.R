#' Estimation of Trends and their Derivatives via Local Polynomial Regression
#'
#' This function is an R function for estimating the trend function
#' and its derivatives in an equidistant time series with local polynomial
#' regression and a fixed bandwidth given beforehand.
#'
#' @param y a numeric vector that contains the time series data ordered from
#' past to present.
#' @param v an integer 0, 1, ... that represents the order of derivative that
#' will be estimated.
#'
#' \tabular{cl}{
#' \strong{Number (v)} \tab \strong{Degree of derivative}\cr
#' \emph{0} \tab The function \emph{f(x)} itself\cr
#' \emph{1} \tab The first derivative \emph{f'(x)}\cr
#' \emph{2} \tab The second derivative \emph{f''(x)}\cr
#' \emph{...} \tab ...
#' }
#' @param p an integer >= (v + 1) that represents the order of polynomial;
#' if p - v is odd, this approach has automatic boundary correction.
#'
#' Examplary for v = 0:
#'
#' \tabular{clcll}{
#' \strong{Number (p)} \tab \strong{Polynomial} \tab
#' \strong{p - v} \tab \strong{p - v odd?} \tab
#' \strong{p usable?}\cr
#' \emph{1} \tab Linear \tab 1 \tab Yes
#' \tab Yes\cr
#' \emph{2} \tab Quadratic \tab 2 \tab No
#' \tab No\cr
#' \emph{3} \tab Cubic \tab 3 \tab Yes
#' \tab Yes\cr
#' \emph{...} \tab ... \tab ... \tab ...
#' \tab ...
#' }
#' @param mu an integer 0, 1, 2, ... that represents the smoothness parameter
#' of the kernel weighting function that will be used; is set to \emph{1} by
#' default.
#'
#' \tabular{cl}{
#'\strong{Number (mu)} \tab \strong{Kernel}\cr
#'\emph{0} \tab Uniform Kernel\cr
#'\emph{1} \tab Epanechnikov Kernel\cr
#'\emph{2} \tab Bisquare Kernel\cr
#'\emph{3} \tab Triweight Kernel\cr
#'\emph{...} \tab ...
#' }
#' @param b a real number 0 < b < 0.5; represents the relative bandwidth that
#' will be used for the smoothing process; is set to \emph{0.15} by default.
#' @param bb can be set to \emph{0} or \emph{1}; the parameter controlling the
#' bandwidth used at the boundary; is set to \emph{0} by default.
#'
#' \tabular{cl}{
#' \strong{Number (bb)} \tab \strong{Estimation procedure at boundary points}\cr
#' \emph{0} \tab Fixed bandwidth on one side with possible large
#' bandwidth on the other side at the boundary\cr
#' \emph{1} \tab The k-nearest neighbor method will be used
#' }
#'
#'@export
#'
#'@details
#'The trend or its derivatives are estimated based on the additive
#'nonparametric regression model for a time series
#'
#'                      y_[t] = m(x_[t]) + eps_[t],
#'
#'where y_[t] is the observed time series, x_[t] is the rescaled time on
#'[0, 1], m(x_[t]) is a smooth trend function and eps_[t] are stationary errors
#'with E(eps_[t]) = 0 (see also Beran and Feng, 2002).
#'
#'This function is part of the package \emph{smoots} and is used in
#'the field of analysing equidistant time series data. It applies the local
#'polynomial regression method to the input data with an arbitrarily
#'selectable bandwidth. By these means, the trend as well as its derivatives
#'can be estimated nonparametrically, even though the result will strongly
#'depend on the bandwidth given beforehand as input.
#'
#'@return The output object is a list with different components:
#'
#'\describe{
#'\item{b}{the chosen (relative) bandwidth; input argument.}
#'\item{bb}{the chosen bandwidth option at the boundaries; input argument.}
#'\item{mu}{the chosen smoothness parameter for the second order kernel; input
#'argument.}
#'\item{n}{the number of observations.}
#'\item{orig}{the original input series; input argument.}
#'\item{p}{the chosen order of polynomial; input argument.}
#'\item{res}{a vector with the estimated residual series; only meaningful for
#'\emph{v = 0}.}
#'\item{v}{the order of derivative; input argument.}
#'\item{ws}{the obtained weighting system matrix for possible advanced
#'analysis.}
#'\item{ye}{a vector with the estimates of the selected nonparametric order of
#'derivative.}
#'}
#'
#'@references
#' Beran, J. and Feng, Y. (2002). Local polynomial fitting with long-memory,
#' short-memory and antipersistent errors. Annals of the Institute of
#' Statistical Mathematics, 54(2), 291-311.
#'
#' Feng, Y., Gries, T. and Fritz, M. (2019). Data-driven
#' local polynomial for the trend and its derivatives in economic time
#' series. Discussion Paper. Paderborn University. (Not yet published)
#'
#' Feng, Y., Gries, T., Letmathe, S. and Schulz, D. (2019). The smoots package
#' in R for semiparametric modeling of trend stationary time series. Discussion
#' Paper. Paderborn University. (Not yet published)
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
#'
#'# Applied gsmooth function for the trend with two different bandwidths
#'results1 <- gsmooth(y, v = 0, p = 1, mu = 1, b = 0.28, bb = 1)
#'results2 <- gsmooth(y, v = 0, p = 1, mu = 1, b = 0.11, bb = 1)
#'trend1 <- results1$ye
#'trend2 <- results2$ye
#'
#'# Plot of the results
#'t <- seq(from = 1947, to = 2019.25, by = 0.25)
#'plot(t, y, type = "l", xlab = "Year", ylab = "log(US-GDP)", bty = "n",
#'     lwd = 2,
#'     main = "Estimated trend for log-quarterly US-GDP, Q1 1947 - Q2 2019")
#'points(t, trend1, type = "l", col = 2, lwd = 1)
#'points(t, trend2, type = "l", col = 4, lwd = 1)
#'legend("bottomright", legend = c("Trend (b = 0.28)", "Trend (b = 0.11)"),
#'      fill = c(2, 4), cex = 0.6)
#'title(sub = expression(italic("Figure 1")), col.sub = "gray47",
#'      cex.sub = 0.6, adj = 0)
#'
#'

# The basic function "gsmooth(y,v,p,mu,b,bb)"-------------------------------

gsmooth <- function(y, v = 0, p = 1, mu = 1, b = 0.15, bb = c(0, 1)) {

  if (((p - v) %% 2) == 0) {
    stop("p - v must be odd.")
  }
  if (all(bb == c(0, 1))) bb <- 0
  n <- length(y)                 # number of observations
  gr <- rep(0, 1, n)             # vector of estimates
  hh <- trunc(n * b + 0.5)       # half bandwidth of the observation time i
  htm <- 2 * hh + 1              # maximal total bandwith of i
  ws <- matrix(0, htm, htm)      # weighting matrix of lpf
  wk <- rep(0, htm)              # kernel weights
  xt <- matrix(0, (p + 1), htm)  # data matrix
  xw <- xt                       # weighted data matrix

  # Linear smoother------------------------------------------------------------

  hr <- c(hh + bb * (hh - 0:hh))  # pointwise right bandwidth of i
  ht <- c(1 + 0:hh) + hr          # pointwise total bandwidth of i

  for(i in (1:(hh + 1))) {

    wk[1:ht[i]] <- (1 - (((1:ht[i]) - i) / (hr[i] + 1)) ** 2) ** (mu)

    for(j in (1:(p + 1))) {

      xt[j, 1:ht[i]] <- (((1:ht[i]) - i) / hh) ** (j - 1)

      # t* = k =: (k - 0.5) / hh is between -1 and 2 for 1 <= k <= ht[i]
	    xw[j,] <- xt[j,] * t(wk)
	  }
    xa <- solve(xw %*% t(xt)) %*% xw
    ws[i,] <- xa[(v+1),]
  }

  # Weighting systems for i > n - hh
  ws[(hh + 2):htm, ] <- (-1) ** (v) * ws[hh:1, htm:1]
  # Normalizing all weighting systems
  ws <- factorial(v) * ws * (n / hh) ** (v)  # normalizing all weighting systems

  ym <- matrix(0,(n - htm + 1), htm)
  for (i in ((hh + 1):(n - hh))) {
    ym[i - hh,] <- y[(i - hh):(i + hh)]  # observations used at point i
  }

    # Estimates for i = 1 -- hh + 1
    gr[1:hh] <- ws[1:hh, ] %*% cbind(y[1:htm])
    # Estimates for hh + 1 < i <= n - hh
    gr[(hh + 1):(n-hh)] <- ym %*% cbind(ws[hh + 1, ])
    # Estimates for i > n - hh
    gr[(n - hh + 1):n] <- ws[(hh + 2):htm,] %*% cbind(y[(n - htm + 1):n])
    result <- list(ye = gr, orig = y, ws = ws, v = v, p = p, mu = mu, b = b,
                   res = y - gr, bb = bb, n = n)

    class(result) = "smoots"  # set package attribute
        attr(result, "function") = "gsmooth"  # set function attribute for
                                                 # own print function

        drop(result)
}
# End of the function
