#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// written by Dominik Schulz; version: 14/04/2020

// Arguments:

// X: The observed time series.
// innov: The residual series that is usually obtained from the 'arima()'
//        function.
// beta: The AR-coefficients ordered from i = 1 to p.
// alpha: The MA-coefficients ordered from j = 1 to q.
// mu: The estimated mean of the series 'X'.
// Fi: The demeaned residual series; this series works as the empirical
//     distribution, from which innovations are generated for the simulation.
// h: The forecasting horizon; this arguments represents the maximum future
//    time period n + h to calculate point forecasts and forecasting intervals
//    for. The calculations are done for all future periods n + 1 to n + h.

// Function for point forecasts with horizon h
// [[Rcpp::export]]
NumericVector fcastCpp(arma::vec X, arma::vec innov,
  arma::vec beta, arma::vec alpha, double mu, int h) {
  int p = beta.size();
  int q = alpha.size();
  NumericVector pq = NumericVector::create(p, q);
  int l = Rcpp::max(pq);
  int n = X.size();
  arma::vec Xsub(l + h);
  arma::vec isub(l + h - 1);
  isub.zeros();
  Xsub.subvec(0, l - 1) = X.subvec(n - l, n - 1);
  isub.subvec(0, l - 1) = innov.subvec(n - l, n - 1);
  arma::rowvec a = (reverse(alpha)).t();
  arma::rowvec b = (reverse(beta)).t();
  for (int i = 0 + l; i < l + h; ++i) {
    Xsub.subvec(i, i) = b * (Xsub.subvec(i - p, i - 1) - mu) +
      a * isub.subvec(i - q, i - 1) + mu;
  }
  arma::vec out = Xsub.subvec(l, l + h - 1);
  return Rcpp::NumericVector(out.begin(), out.end());
}

// Function for simulated 'true' forecasts
// [[Rcpp::export]]
NumericVector tfcastCpp(arma::vec X, arma::vec innov, NumericVector Fi,
  arma::vec beta, arma::vec alpha, double mu, int h) {
  int p = beta.size();
  int q = alpha.size();
  NumericVector pq = NumericVector::create(p, q);
  int l = Rcpp::max(pq);
  int n = X.size();
  arma::vec Xsub(l + h);
  arma::vec isub(l + h );
  arma::vec sampleInnov = Rcpp::sample(Fi, h, true);
  Xsub.subvec(0, l - 1) = X.subvec(n - l, n - 1);
  isub.subvec(0, l - 1) = innov.subvec(n - l, n - 1);
  isub.subvec(l, l + h - 1) = sampleInnov;
  arma::rowvec a = (reverse(alpha)).t();
  arma::rowvec b = (reverse(beta)).t();
  for (int i = 0 + l; i < l + h; ++i) {
    Xsub.subvec(i, i) = b * (Xsub.subvec(i - p, i - 1) - mu) +
      a * isub.subvec(i - q, i - 1) + mu + isub(i);
  }
  arma::vec out = Xsub.subvec(l, l + h - 1);
  return Rcpp::NumericVector(out.begin(), out.end());
}

// Function for ARMA simulation via bootstrap
// [[Rcpp::export]]
NumericVector arimaSimBootCpp(NumericVector Fi,
  arma::vec beta, arma::vec alpha, double mu, int nStart) {
  int p = beta.size();
  int q = alpha.size();
  NumericVector pq = NumericVector::create(p, q);
  int l = Rcpp::max(pq);
  int n = Fi.size();
  arma::vec Xsub(l + n + nStart);
  Xsub.zeros();
  arma::vec sampleEps = Rcpp::sample(Fi, n + nStart, true);
  arma::vec isub = Xsub;
  isub.subvec(l, l + nStart + n - 1) = sampleEps;
  arma::rowvec a = (reverse(alpha)).t();
  arma::rowvec b = (reverse(beta)).t();
  for (int i = 0 + l; i < l + n + nStart; ++i) {
    Xsub.subvec(i, i) = b * Xsub.subvec(i - p, i - 1) + a * isub.subvec(i - q, i - 1) +
      isub(i);
  }
  arma::vec out = Xsub.subvec(l + nStart, l + nStart + n - 1) + mu;
  return Rcpp::NumericVector(out.begin(), out.end());
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R

*/
