#' Estimation of Nonparametric Trend Functions via Kernel Regression
#'
#' This function estimates the nonparametric trend function in an equidistant
#' time series with Nadaraya-Watson kernel regression.
#'
#' @param y a numeric vector that contains the time series data ordered from past to present.
#' @param mu an integer 0, 1, 2, ... that represents the smoothness parameter
#' of the second order kernel function that will be used; is set to \emph{1} by
#' default.
#'
#'\tabular{cl}{
#'\strong{Number (mu)} \tab \strong{Kernel}\cr
#'\emph{0} \tab Uniform Kernel\cr
#'\emph{1} \tab Epanechnikov Kernel\cr
#'\emph{2} \tab Bisquare Kernel\cr
#'\emph{3} \tab Triweight Kernel\cr
#'\emph{...} \tab ...
#' }
#' @param b a real number 0 < b < 0.5; represents the relative bandwidth that
#' will be used for the smoothing process; is set to \emph{0.15} by default.
#' @param bb can be set to 0 or 1; the parameter controlling the bandwidth used
#' at the boundary; is set to \emph{0} by default.
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
#'
#'The trend or its derivatives are estimated based on the additive
#'nonparametric regression model for a time series
#'
#'                      y_[t] = m(x_[t]) + eps_[t],
#'
#'where y_[t] is the observed time series, x_[t] is the rescaled time on
#'[0, 1], m(x_[t]) is a smooth trend function and eps_[t] are stationary errors
#'with E(eps_[t]) = 0.
#'
#'This function is part of the package \emph{smoots} and is used for
#'the estimation of trends in equidistant time series. The applied method
#'is a kernel regression with arbitrarily selectable second order
#'kernel, relative bandwidth and boundary method. Especially the chosen
#'bandwidth has a strong impact on the final result and has thus to be
#'selected carefully. This approach is not recommended by the authors of this
#'package.
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
#'\item{res}{a vector with the estimated residual series.}
#'\item{ye}{a vector with the estimates of the nonparametric trend.}
#'}
#'
#'@references
#'Feng, Y. (2009). Kernel and Locally Weighted Regression. Verlag für
#'Wissenschaft und Forschung, Berlin.
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
#'#Applied knmooth function for the trend with two different bandwidths
#'trend1 <- knsmooth(y, mu = 1, b = 0.28, bb = 1)$ye
#'trend2 <- knsmooth(y, mu = 1, b = 0.05, bb = 1)$ye
#'
#'# Plot of the results
#'t <- seq(from = 1947, to = 2019.25, by = 0.25)
#'plot(t, y, type = "l", xlab = "Year", ylab = "log(US-GDP)", bty = "n",
#'     lwd = 2,
#'     main = "Estimated trend for log-quarterly US-GDP, Q1 1947 - Q2 2019")
#'points(t, trend1, type = "l", col = 2, lwd = 1)
#'points(t, trend2, type = "l", col = 4, lwd = 1)
#'legend("bottomright", legend = c("Trend (b = 0.28)", "Trend (b = 0.05)"),
#'      fill = c(2, 4), cex = 0.6)
#'title(sub = expression(italic("Figure 1")), col.sub = "gray47",
#'      cex.sub = 0.6, adj = 0)
#'
#'

knsmooth <- function(y, mu = 1, b = 0.15, bb = c(0, 1)) {

  if (all(bb == c(0, 1))) bb = 0
  n <- length(y)           #number of observations
  hn <- trunc(b * n + 0.5)
  gr <- rep(0, n)         #vector of estimates
  if (hn > trunc((n - 1) / 2)) {hn <- trunc((n - 1) / 2)} ### maximal bandwidth

  ####### kernel smoothing
  hr <- c(hn + bb * (hn - (0:hn)))    #pointwise right bandwith of i
  ht <- c(1 + (0:hn)) + hr       #pointwise total bandwith of i

  for (i in (1:(hn + 1))) {
    u <- ((1:ht[i]) - i) / (hr[i] + 0.5)

    ###### The kernel functions
    wk <- 1 / 2 * (1 - u^2)^mu
    wk <- wk / sum(wk) ##### so that the sum of all weights is 1.

    if (i <= hn) {
      gr[i] <- sum(wk * y[1:ht[i]])
      gr[n - i + 1] <- sum(wk * y[n:(n - ht[i] + 1)])
    }

    if(i == hn + 1) {

      for(j in i:(n-hn)) {
        gr[j] <- sum(wk * y[(j - hn):(j + hn)])
      }
    }
  }

  result <- list(ye = gr, orig = y, res = y - gr, mu = mu, b = b, bb = bb,
                 n = n)

  class(result) = "smoots"  # set package attribute
  attr(result, "function") = "knsmooth"  # set function attribute for
  # own print function

  return(result)  ####The results are given in gr
}
