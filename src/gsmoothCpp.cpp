#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// FUnctions written by Dominik Schulz, 12/05/2020

// Function to help obtain a sequence from 'from' to 'to' by step 1
// in C++
// [[Rcpp::export]]
arma::vec seqCpp(int from, int to){
  int n = to - from + 1;
  arma::vec seqOut(n);
  for (int i = from; i < to + 1; ++i) {
    seqOut(i - from) = i;
  }
  return seqOut;
}

// Function to help obtain a sequence from 'from' to 'to' by step 1
// in C++ (returns a rowvec instead)
// [[Rcpp::export]]
arma::rowvec rseqCpp(int from, int to) {
  int n = to - from + 1;
  arma::rowvec seqOut(n);
  for (int i = from; i < to + 1; ++ i) {
    seqOut(i - from) = i;
  }
  return seqOut;
}

// Function to calculate the factorial of an 'int' in C++
// [[Rcpp::export]]
int factorialCpp(int k) {
  int fac = 1;
  if (k > 1) {
    for (int i = 2; i < k + 1; ++i) {
      fac *= i;
    }
  }
  return fac;
}

// C++ version of gsmoothCalc
// [[Rcpp::export]]
NumericVector gsmoothCalcCpp(arma::vec y, int v, int p, int mu, double b, int bb) {
  int n = y.size();
  arma::vec gr(n);
  int hh = std::trunc(n * b + 0.5);
  int htm = 2 * hh + 1;
  arma::mat ws(htm, htm);
  arma::rowvec wk(htm);
  wk.zeros();
  arma::mat xt = arma::zeros((p + 1), htm);
  arma::mat xw = xt;

  arma::vec seqhh = seqCpp(0, hh);
  arma::vec hr = hh + bb * (hh - seqhh);
  arma::vec ht = 1 + seqhh + hr;
  arma::mat xa = xt;

  arma::rowvec sequ;

  arma::rowvec sequhh;

  for (int i = 0; i < hh + 1; ++i) {

    sequ = rseqCpp(1, ht(i)) - (i + 1);

    wk.subvec(0, ht(i) - 1) = arma::pow(1.0 - arma::pow(sequ / (hr(i) + 1), 2.0), mu);

    sequhh = sequ / hh;

    for (int j = 0; j < p + 1; ++j) {
      xt.submat(j, 0, j, ht(i) - 1) = arma::pow(sequhh, j);
      xw.row(j) = xt.row(j) % wk;
    }

    xa.submat(0, 0, p, htm - 1) = arma::solve(xw * xt.t(), xw);

    ws.row(i) = xa.row(v);
  }

  ws.rows(hh + 1, htm - 1) = std::pow(-1.0, v) * arma::flipud(arma::fliplr(ws.submat(0, 0, hh - 1, htm - 1)));
  ws = factorialCpp(v) * ws * std::pow(n / double(hh), v);
  arma::mat ym(n - htm + 1, htm);
  for (int i = hh; i < n - hh; ++i) {
    ym.row(i - hh) = (y.subvec(i - hh, i + hh)).t();
  }
  gr.subvec(0, hh - 1) = ws.rows(0, hh - 1) * y.subvec(0, htm - 1);
  gr.subvec(hh, n - hh - 1) = ym * (ws.row(hh)).t();
  gr.subvec(n - hh, n - 1) = ws.rows(hh + 1, htm - 1) * y.subvec(n - htm, n - 1);

  return Rcpp::NumericVector(gr.begin(), gr.end());
}

// C++ version of gsmoothCalc2
// [[Rcpp::export]]
List gsmoothCalc2Cpp(arma::vec y, int v, int p, int mu, double b, int bb) {
  int n = y.size();
  arma::vec gr(n);
  int hh = std::trunc(n * b + 0.5);
  int htm = 2 * hh + 1;
  arma::mat ws(htm, htm);
  arma::rowvec wk(htm);
  wk.zeros();
  arma::mat xt = arma::zeros((p + 1), htm);
  arma::mat xw = xt;

  arma::vec seqhh = seqCpp(0, hh);
  arma::vec hr = hh + bb * (hh - seqhh);
  arma::vec ht = 1 + seqhh + hr;
  arma::mat xa = xt;

  arma::rowvec sequ;
  arma::rowvec sequhh;

  for (int i = 0; i < hh + 1; ++i) {

    sequ = rseqCpp(1, ht(i)) - (i + 1);
    wk.subvec(0, ht(i) - 1) = arma::pow(1.0 - arma::pow(sequ / (hr(i) + 1), 2.0), mu);
    sequhh = sequ / hh;

    for (int j = 0; j < p + 1; ++j) {
      xt.submat(j, 0, j, ht(i) - 1) = arma::pow(sequhh, j);
      xw.row(j) = xt.row(j) % wk;
    }

    xa.submat(0, 0, p, htm - 1) = arma::solve(xw * xt.t(), xw);

    ws.row(i) = xa.row(v);
  }

  ws.rows(hh + 1, htm - 1) = std::pow(-1.0, v) * arma::flipud(arma::fliplr(ws.submat(0, 0, hh - 1, htm - 1)));
  ws = factorialCpp(v) * ws * std::pow(n / double(hh), v);
  arma::mat ym(n - htm + 1, htm);
  for (int i = hh; i < n - hh; ++i) {
    ym.row(i - hh) = (y.subvec(i - hh, i + hh)).t();
  }
  gr.subvec(0, hh - 1) = ws.rows(0, hh - 1) * y.subvec(0, htm - 1);
  gr.subvec(hh, n - hh - 1) = ym * (ws.row(hh)).t();
  gr.subvec(n - hh, n - 1) = ws.rows(hh + 1, htm - 1) * y.subvec(n - htm, n - 1);

  List listOut = List::create(_["ye"] = Rcpp::NumericVector(gr.begin(), gr.end()),
                              _["ws"] = ws);

  return listOut;

}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R

*/
