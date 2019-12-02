#' Advanced Data-driven Nonparametric Regression for the Trend in Equidistant
#' Time Series
#'
#' This function runs an iterative plug-in algorithm to find the optimal
#' bandwidth for the estimation of the nonparametric trend in equidistant
#' time series (with short-memory errors) and then employs the resulting
#' bandwidth via either local polynomial or kernel regression. This function
#' allows for more flexibility in its arguments than \emph{msmooth}.
#'
#'@param y a numeric vector that contains the time series ordered from past to present.
#'@param p an integer 1 (local linear regression) or 3 (local cubic regression);
#'represents the order of polynomial within the local polynomial regression
#'(see also the 'Details' section); is set to \emph{1} by default;
#'is automatically set to \emph{1} if \emph{method = "kr"}.
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
#'@param Mcf  method for estimating the variance factor c_[f] by the IPI (see
#'also the 'Details' section); is set to \emph{NP} by default. \cr
#'\tabular{cl}{
#'\strong{Method} \tab \strong{Explanation}\cr
#'\emph{NP} \tab Nonparametric estimation using the Bartlett window \cr
#'\emph{ARMA} \tab Estimation on the assumption that the residuals follow an
#'ARMA model\cr
#'\emph{AR} \tab Estimation on the assumption that the residuals follow an
#'AR model\cr
#'\emph{MA} \tab Estimation on the assumption that the residuals follow an
#'MA model
#'}
#'@param InfR a character object that represents the inflation
#'rate in the form h_[d] = h^[a] for the bandwidth in the estimation of
#'I[m^(k)] (see also the 'Details' section); is set to \emph{Opt} by default.
#'
#'\tabular{cl}{
#'\strong{Inflation rate} \tab \strong{Description}\cr
#'\emph{Opt} \tab Optimal inflation rate a_[p,O] (5/7 for p = 1; 9/11 for
#'p = 3)\cr
#'\emph{Nai} \tab Naive inflation rate a_[p,N] (5/9 for p = 1; 9/13 for p = 3)\cr
#'\emph{Var} \tab Stable inflation rate a_[p,S] (1/2 for p = 1 and p = 3)
#'}
#'@param bStart a numeric object that indicates the starting value of the
#'bandwidth for the iterative process; should be 0 < bStart < 0.5; is set to
#'\emph{0.15} by default.
#'@param bvc a character object that indicates whether an enlarged bandwidth is
#'being used for the estimation of the variance factor c_[f] (see also the
#''Details' section) or not; is set to \emph{Y} by default.
#'
#'\tabular{cl}{
#'\strong{Bandwidth} \tab \strong{Description}\cr
#'\emph{Y} \tab Using an enlarged bandwidth\cr
#'\emph{N} \tab Using a bandwidth without enlargement
#'}
#'@param bb can be set to \emph{0} or \emph{1}; the parameter controlling the
#'bandwidth used at the boundary; is set to \emph{1} by default.
#'
#'\tabular{cl}{
#'\strong{Number (bb)} \tab \strong{Estimation procedure at boundary points}\cr
#'\emph{0} \tab Fixed bandwidth on one side with possible large
#'bandwidth on the other side at the boundary\cr
#'\emph{1} \tab The k-nearest neighbor method will be used
#'}
#'@param cb a numeric value that indicates the percentage of omitted
#'observations on each side of the observation period for the automated
#'bandwidth selection; is set to \emph{0.05} by default.
#'@param method the final smoothing approach; \emph{"lpr"} represents the local
#'polynomial regression, whereas \emph{"kr"} implements a kernel regression;
#'is set to \emph{"lpr"} by default.
#'
#'@details
#'
#'A nonparametric regression of the trend in a time series is based on the
#'additive model
#'
#'                      y_[t] = m(x_[t]) + eps_[t],
#'
#'where y_[t] is the observed time series in question, x_[t] is the rescaled
#'time on [0, 1], m(x_[t]) is the nonparametric trend and eps_[t] are the
#'errors with E(eps_[t]) = 0 (see also Beran and Feng, 2002). With this
#'function m(x_[t]) can be estimated without a parametric model assumption for
#'the error series. Thus, after estimating and removing the trend, any suitable
#'parametric model, e.g. an ARMA(p, q) model, can be fitted to the residuals.
#'
#'The iterative-plug-in (IPI) algorithm, which numerically minimizes the
#'Asymptotic Mean Squared Error (AMISE), was proposed by Feng, Gries and Fritz
#'(2019).
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
#'The function calculates suitable estimates for c_[f], the variance factor,
#'and I[m^(k)] over different iterations. In each iteration, a bandwidth is
#'obtained in accordance with the AMISE that once more serves as an input for
#'the following iteration. The process repeats until either convergence or the
#'40th iteration is reached. For further details on the asymptotic theory or
#'the algorithm, please consult Feng, Gries and Fritz (2019) or Feng et al.
#'(2019).
#'
#'To apply the function, more arguments are needed compared to the similar
#'function \emph{msmooth}: a data input \emph{y}, an order of polynomial
#'\emph{p}, a kernel weighting function defined by the smoothness parameter
#'\emph{mu}, a variance factor estimation method \emph{Mcf}, an inflation rate
#'setting \emph{InfR} (see also Beran and Feng, 2002), a starting value for the
#'relative bandwidth \emph{bStart}, an inflation setting for the variance
#'factor estimation \emph{bvc}, a boundary method \emph{bb}, a boundary cut-off
#'percentage \emph{cb} and a final smoothing method \emph{method}.
#'In fact, aside from the input vector \emph{y}, every argument has a default
#'setting that can be adjusted for the individual case. Theoretically, the
#'initial bandwidth does not affect the selected optimal bandwidth. However, in
#'practice local minima of the AMISE might exist and influence the selected
#'bandwidth. Therefore, the default setting is \emph{bStart = 0.15}, which is a
#'compromise between the starting values \emph{bStart = 0.1} for \emph{p = 1}
#'and \emph{bStart = 0.2} for \emph{p = 3} that were proposed by Feng, Gries
#'and Fritz (2019). In the rare case of a clearly unsuitable optimal bandwidth,
#'a starting bandwidth that differs from the default value is a first
#'possible approach to obtain a better result. Other argument adjustments can
#'be tried as well. For more specific information on the input arguments
#'consult the section \emph{Arguments}.
#'
#'When applying the function, an optimal bandwidth is obtained based on the
#'IPI algorithm proposed by Feng, Gries and Fritz (2019). In a second step,
#'the nonparametric trend of the series is calulated with respect
#'to the chosen bandwidth and the selected regression method (\emph{lpf} or
#'\emph{kr}). Please note that \emph{method = "lpf"} is strongly recommended by
#'the authors. Moreover, it is notable that \emph{p} is automatically set to 1
#'for \emph{method = "kr"}. The output object is then a list that contains,
#'among other components, the original time series, the estimated trend values
#'and the series without the trend.
#'
#'The default print method for this function delivers only key numbers such as
#'the iteration steps and the generated optimal bandwidth rounded to the fourth
#'decimal. The exact numbers and results such as the estimated nonparametric
#'trend series are saved within the output object and can be addressed via the
#'\emph{$} sign.
#'
#'@return The function returns a list with different components:
#'
#'\describe{
#'\item{AR.BIC}{the Bayesian Information Criterion of the optimal AR(p) model
#'when estimating the variance factor via autoregressive models (if
#'calculated; calculated for \emph{Mcf = "AR"}).}
#'\item{ARMA.BIC}{the Bayesian Information Criterion of the optimal ARMA(p,q)
#'model when estimating the variance factor via autoregressive-moving-average
#'models (if calculated; calculated for \emph{Mcf = "ARMA"}).}
#'\item{cb}{the percentage of omitted observations on each side of the
#'observation period; input argument.}
#'\item{b0}{the optimal bandwidth chosen by the IPI-algorithm.}
#'\item{bb}{the boundary bandwidth method used within the IPI; input argument.}
#'\item{bStart}{the starting value of the (relative) bandwidth; input
#'argument.}
#'\item{bvc}{indicates whether an enlarged bandwidth was used for the variance
#'factor estimation or not; input argument.}
#'\item{cf0}{the estimated variance factor.}
#'\item{cf0.AR}{the estimated variance factor obtained by estimation of
#'autoregressive models (if calculated; calculated for
#'\emph{Mcf = "AR"}).}
#'\item{cf0.ARMA}{the estimated variance factor obtained by estimation of
#'autoregressive-moving-average models (if calculated; calculated for
#'\emph{Mcf = "ARMA"}).}
#'\item{cf0.LW}{the estimated variance factor obtained by Lag-Window Spectral
#'Density Estimation following Bühlmann (1996) (if calculated; calculated for
#'\emph{Mcf = "NP"}).}
#'\item{cf0.MA}{the estimated variance factor obtained by estimation of
#'moving-average models (if calculated; calculated for
#'\emph{Mcf = "MA"}).}
#'\item{I2}{the estimated value of I[m(k)].}
#'\item{InfR}{the setting for the inflation rate; input argument.}
#'\item{iterations}{the bandwidths of the single iterations steps.}
#'\item{L0.opt}{the optimal bandwidth for the lag-window spectral density
#'estimation (if calculated).}
#'\item{MA.BIC}{the Bayesian Information Criterion of the optimal MA(q) model
#'when estimating the variance factor via moving-average models (if
#'calculated; calculated for \emph{Mcf = "MA"}).}
#'\item{Mcf}{the estimation method for the variance factor estimation; input
#'argument.}
#'\item{mu}{the smoothness parameter of the second order kernel; input
#'argument.}
#'\item{n}{the number of observations.}
#'\item{niterations}{the total number of iterations until convergence.}
#'\item{orig}{the original input series; input argument.}
#'\item{p.BIC}{the order p of the optimal AR(p) or ARMA(p,q) model when
#'estimating the variance factor via autoregressive or autoregressive-moving
#'average models (if calculated; calculated for \emph{Mcf = "AR"} and
#'\emph{Mcf = "ARMA"}).}
#'\item{p}{the order of polynomial used in the IPI-algorithm; also used for the
#'final smoothing, if \emph{method = "lpr"}; input argument.}
#'\item{q.BIC}{the order q of the optimal MA(q) or ARMA(p,q) model when
#'estimating the variance factor via moving-average or autoregressive-moving
#'average models (if calculated; calculated for \emph{Mcf = "MA"} and
#'\emph{Mcf = "ARMA"}).}
#'\item{res}{the estimated residual series.}
#'\item{ws}{the weighting systems used within the local polynomial regression;
#'only exists, if the final smoothing is done via a local polynomial
#'regression.}
#'\item{ye}{the nonparametric estimates of the trend.}
#'}
#'
#'@export
#'
#'@references
#' Beran, J. and Feng, Y. (2002). Local polynomial fitting with long-memory,
#' short-memory and antipersistent errors. Annals of the Institute of
#' Statistical Mathematics, 54(2), 291-311.
#'
#' Bühlmann, P. (1996). Locally adaptive lag-window spectral estimation.
#' Journal of Time Series Analysis, 17(3), 247-270.
#'
#' Feng, Y., Gries, T. and Fritz, M. (2019). Data-driven
#' local polynomial for the trend and its derivatives in economic time
#' series. Discussion Paper. Paderborn University. (Not yet pubslished)
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
#'### Example 1: US-GDP ###
#'
#'# Logarithm of test data
#'test_data <- gdpUS
#'y <- log(test_data$GDP)
#'
#'# Applied tsmooth function for the trend
#'result <- tsmooth(y, p = 1, mu = 1, Mcf = "NP", InfR = "Opt",
#'                      bStart = 0.1, bvc = "Y")
#'trend1 <- result$ye
#'
#'# Plot of the results
#'t <- seq(from = 1947, to = 2019.25, by = 0.25)
#'plot(t, y, type = "l", xlab = "Year", ylab = "log(US-GDP)", bty = "n",
#'     lwd = 1, lty = 3,
#'     main = "Estimated trend for log-quarterly US-GDP, Q1 1947 - Q2 2019")
#'points(t, trend1, type = "l", col = 2, lwd = 1)
#'title(sub = expression(italic("Figure 1")), col.sub = "gray47",
#'      cex.sub = 0.6, adj = 0)
#'result
#'
#'\donttest{
#'### Example 2: German Stock Index ###
#'
#'# Obtain log-transformation of the returns
#'returns <- diff(log(dax$Close))
#'rt <- returns - mean(returns)
#'
#'# Apply 'smoots' function to the log-data
#'# -> the log-transformation is assumed to follow the additive model
#'yt <- log(rt^2)
#'
#'# In this case, the optimal inflation rate is used for p = 3.
#' est <- tsmooth(yt, p = 3, InfR = "Opt")
#' m_xt <- est$ye
#'
#'# Obtain the standardized returns 'eps' and the scale function 's'
#' sqrtC_s <- exp(m_xt / 2)
#' eps_sqrtC <- rt / sqrtC_s
#' C <- 1 / var(eps_sqrtC)
#' eps <- eps_sqrtC * sqrt(C)
#' s <- sqrtC_s / sqrt(C)
#'
#'# -> 'eps' can now be analyzed by any suitable GARCH-type model.
#'#    The total volatilities are then the product of the conditional
#'#    volatilities obtained from 'eps' and the scale function 's'.
#'}

# The main function------------------------------------------------------------

tsmooth <- function(y, p = c(1, 3), mu = c(0, 1, 2, 3),
                    Mcf = c("NP", "ARMA", "AR", "MA"),
                    InfR = c("Opt", "Nai", "Var"), bStart = 0.15,
                    bvc = c("Y", "N"), bb = c(0, 1), cb = 0.05,
                    method = c("lpr", "kr")) {

     # Input if no inputs were made to specific arguments
     if (all(mu == c(0, 1, 2, 3))) mu <- 1
     if (all(Mcf == c("NP", "ARMA", "AR", "MA"))) Mcf <- "NP"
     if (all(InfR == c("Opt", "Nai", "Var"))) InfR <- "Opt"
     if (all(bvc == c("Y", "N"))) bvc <- "Y"
     if (all(bb == c(0, 1))) bb <- 1
     if (all(method == c("lpr", "kr"))) method <- "lpr"
     if (all(p == c(1, 3)) || method == "kr") p <- 1

     # Stop, if incorrect inputs were given
     if (!(p %in% c(1, 3))) {
       stop("Input of argument 'p' incorrect. It must be either 1 or 3.")
     }
     if (!(mu %in% 0:4)) {
       stop("Input of argument 'mu' incorrect. It must be an integer from 0 to 4.")
     }
     if (!(Mcf %in% c("NP", "ARMA", "AR", "MA"))) {
       stop("Input of argument 'Mcf' incorrect. Method not recognized.")
     }
     if (!(InfR %in% c("Opt", "Nai", "Var"))) {
       stop("Input of argument 'InfR' incorrect. Input not recognized.")
     }
     if (!(InfR %in% c("Opt", "Nai", "Var"))) {
       stop("Input of argument 'InfR' incorrect. Input not recognized.")
     }
     if (!(bvc %in% c("Y", "N"))) {
       stop("Input of argument 'bvc' incorrect. Input not recognized.")
     }
     if(!(is.numeric(bb) && bb %in% c(0, 1))) {
       stop("Input of argument 'bb' incorrect. It must be either 0 or 1.")
     }
     if(!(is.numeric(cb) && cb >= 0 && cb < 0.5)) {
       stop("Input of argument 'cb' incorrect. It must be a numeric value
            between 0 and 0.5.")
     }
     if (!(method %in% c("lpr", "kr"))) {
       stop("Input of argument 'method' incorrect. Method not recognized.")
     }
     if (!(is.numeric(bStart) && length(bStart) == 1)) {
       stop("Input of argument 'bStart' incorrect.
            One numerical value between 0 and 0.5 is needed.")
     }
     if (bStart <= 0 || bStart >= 0.5) {
       message("NOTE: 'bStart' between 0 and 0.5 is recommended.")
     }

     # using the knn idea, bb=1, or not
     # cb=x, 0, 0.025 or 0.05, x*n obs at each end not be used for choosing b
     # In this code ccf and bStart will be provided by the user
     # Input parameters
     n <- length(y)
     k <- p + 1
     pd <- p + 2
     runc <- 1
     n1 <- trunc(n * cb)

     # New method for the kernel constants with p = 1 or 3
     # Kernel constants
     m <- 1000000  # for the numerical integral
     u <- (-m:m) / (m + 0.5)
     # For p=1, any kernel in the given c-mu class is possible
     if (p == 1) {
       wkp <- (1 - u^2)^(mu)  # a standardization is not necessary
     }

     # For p=3, the four 4-th order kernels (Table 5.7, Mueller, 1988)
     # table saved internally
     if (p == 3) {  # the constant factor does not play any role
       wkp <- lookup$p3_lookup[mu + 1][[1]](u)
     }
                    Rp <- sum(wkp^2) / m
                    mukp <- sum((u^k) * wkp) / m

     # Two constant in the bandwidth
     c1 <- factorial(k)^2 / (2 * k)
     c2 <- (1 - 2 * cb) * Rp / (mukp)^2

     steps <- rep(NA, 40)
     bd_func <- lookup$InfR_lookup[as.character(p), InfR][[1]]
     bv_func <- lookup$bvc_lookup[as.character(p), mu + 1][[1]]

     # The main iteration------------------------------------------------------

     noi <- 40  # maximal number of iterations

     for (i in 1:noi) {
       if (runc == 1) {
         if (i > 1) {bold1 <- bold}
         if (i == 1) {bold <- bStart} else {bold <- bopt}

         # Look up the EIM inflation rate in the internal list 'lookup'
         bd <- bd_func(bold)

         if (bd >= 0.49) {bd <- 0.49}
         yed <- gsmooth(y, k, pd, mu, bd, bb)$ye
         I2 <- sum(yed[(n1 + 1):(n - n1)]^2) / (n - 2 * n1)

         # Use an enlarged bandwidth for estimating cf or not (Feng/Heiler, 2009)
         if (bvc == "Y") {

           # Look up the enlarged bandwidth
           bv <- bv_func(bold)

         } else {bv <- bold}

         if (bv >= 0.49) {bv <- 0.49}

         ye <- gsmooth(y, 0, p, mu, bv, bb)$ye

         # The optimal bandwidth-----------------------------------------------

         # Estimating the variance factor
         yd <- y - ye

         # Lag-Window for estimating the variance
         if (Mcf == "NP") {
           cf0.LW.est <- cf0.LW.est(yd)  #[(n1+1):(n-n1)])
           cf0.LW <- cf0.LW.est$cf0.LW
           cf0 <- cf0.LW
           cf0.AR <- NA
           cf0.MA <- NA
           cf0.ARMA <- NA
           L0.opt <- cf0.LW.est$L0.opt
           ARMA.BIC <- NA
           AR.BIC <- NA
           MA.BIC <- NA
           p.BIC <- NA
           q.BIC <- NA
         } else if (Mcf == "ARMA") {
           cf0.ARMA.est <- cf0.ARMA.est(yd)
           cf0.ARMA <- cf0.ARMA.est$cf0.ARMA
           cf0.AR <- NA
           cf0.MA <- NA
           cf0 <- cf0.ARMA
           AR.BIC <- NA
           MA.BIC <- NA
           ARMA.BIC <- cf0.ARMA.est$ARMA.BIC
           p.BIC <- cf0.ARMA.est$p.BIC
           q.BIC <- cf0.ARMA.est$q.BIC
           cf0.LW <- NA
           L0.opt <- NA
         } else if (Mcf == "AR") {
           cf0.AR.est <- cf0.AR.est(yd)  # [(n1+1):(n-n1)])
           cf0.AR <- cf0.AR.est$cf0.AR
           cf0 <- cf0.AR
           cf0.LW <- NA
           cf0.MA <- NA
           cf0.ARMA <- NA
           AR.BIC <- cf0.AR.est$AR.BIC
           p.BIC <- cf0.AR.est$p.BIC
           ARMA.BIC <- NA
           MA.BIC <- NA
           q.BIC <- NA
           L0.opt <- NA
         } else if (Mcf == "MA") {
           cf0.MA.est <- cf0.MA.est(yd)
           cf0.MA <- cf0.MA.est$cf0.MA
           cf0 <- cf0.MA
           cf0.LW <- NA
           cf0.AR <- NA
           cf0.ARMA <- NA
           MA.BIC <- cf0.MA.est$MA.BIC
           q.BIC <- cf0.MA.est$q.BIC
           ARMA.BIC <- NA
           AR.BIC <- NA
           p.BIC <- NA
           L0.opt <- NA
         }

	       c3 = cf0 / I2

         if (p == 1) {
           bopt = (c1 * c2 * c3)^(1 / 5) * n^(-1 / 5)
           if (bopt < n^(-5 / 7)) {bopt = n^(-5 / 7)}
         }
         if (p == 3) {
           bopt = (c1 * c2 * c3)^(1 / 9) * n^(-1 / 9)
         if (bopt < n^(-9 / 11)) {bopt = n^(-9 / 11)}
         }
           if (bopt > 0.49) {bopt = 0.49}
	         steps[i] = bopt
                 if (i > 2 && abs(bold - bopt) / bopt < 1 / n) {runc = 0}
                 if (i > 3 && abs(bold1 - bopt) / bopt < 1 / n) {
                               bopt = (bold + bopt) / 2
                               runc = 0
                 }
       }
     }

     # Smooth with the selected bandwidth--------------------------------------

     if (p == 1 && bopt < n^(-5 / 7)) {bopt <- n^(-5 / 7)}
     if (p == 3 && bopt < n^(-9 / 11)) {bopt <- n^(-9 / 11)}
     if (bopt > 0.49) {bopt <- 0.49}

     if (method == "lpr") {
       est.opt <- gsmooth(y, 0, p, mu, bopt, bb)
       ye <- est.opt$ye
       ws <- est.opt$ws
       res <- y - ye
       results <- list(ye = est.opt$ye, b0 = bopt, cf0 = cf0, I2 = I2,
                       cf0.AR = cf0.AR, cf0.MA = cf0.MA, cf0.ARMA = cf0.ARMA,
                       cf0.LW = cf0.LW, L0.opt = L0.opt,
                       AR.BIC = AR.BIC, MA.BIC = MA.BIC, ARMA.BIC = ARMA.BIC,
                       p.BIC = p.BIC, q.BIC = q.BIC, n = n,
                       niterations = length(steps[!is.na(steps)]), res = res,
                       orig = y, iterations = steps[!is.na(steps)],
                       p = p, mu = mu, Mcf = Mcf, InfR = InfR, bStart = bStart,
                       bvc = bvc, bb = bb, cb = cb)
       attr(results, "method") = "lpr"
     } else if (method == "kr") {
       est.opt <- knsmooth(y, mu, bopt, bb)
       res <- est.opt$res
       results <- list(ye = est.opt$ye, b0 = bopt, cf0 = cf0, I2 = I2,
                       cf0.AR = cf0.AR, cf0.MA = cf0.MA, cf0.ARMA = cf0.ARMA,
                       cf0.LW = cf0.LW, L0.opt = L0.opt,
                       AR.BIC = AR.BIC, MA.BIC = MA.BIC, ARMA.BIC = ARMA.BIC,
                       p.BIC = p.BIC, q.BIC = q.BIC, n = n,
                       niterations = length(steps[!is.na(steps)]), res = res,
                       orig = y, iterations = steps[!is.na(steps)],
                       p = p, mu = mu, Mcf = Mcf, InfR = InfR, bStart = bStart,
                       bvc = bvc, bb = bb, cb = cb)
       attr(results, "method") = "kr"
     }

         class(results) <- "smoots"
         attr(results, "function") = "tsmooth"
         results
}
# End of the function
