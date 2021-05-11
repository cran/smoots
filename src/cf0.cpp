#include <Rcpp.h>
using namespace Rcpp;

// Autocovariance calculation from given data
// [[Rcpp::export]]
NumericVector acovCpp(NumericVector Xt, int lagMax) {
  double meanX = Rcpp::mean(Xt);
  NumericVector acovOut(lagMax + 1);
  int n = Xt.size();
  for (int i = 0; i < lagMax + 1; ++i) {
    acovOut[i] = Rcpp::sum((Xt[Range(0, n - 1 - i)] - meanX) * (Xt[Range(i, n - 1)] - meanX)) / n;
  }
  return acovOut;
}

// C++ equivalent for seq() with NumericVector as data type
// [[Rcpp::export]]
NumericVector seqNum(int begin, int end) {
  IntegerVector series = Rcpp::seq(begin, end);
  NumericVector out = as<NumericVector>(series);
  return out;
}

// Buehlmann algorithm for estimating cf0
// [[Rcpp::export]]
List cf0Cpp(NumericVector Xt) {
  int n = Xt.size();
  NumericVector ga = acovCpp(Xt, n - 1);
  int nit = 20;
  int runc = 1;
  NumericVector L(nit + 1);
  L[0] = std::trunc(n / 2.0 + 0.5);
  double c1 = (std::pow(ga[0], 2.0) + 2.0 * Rcpp::sum(Rcpp::pow(ga[Range(1, n - 1)], 2.0))) / (4.0 * M_PI);
  double c2 = 0;
  double c10 = 0;
  double c20 = 0;
  int L1 = 0;
  int LGopt = 0;

  for(int i = 0; i < nit; i++) {
    if(runc == 1) {
      L1 = std::trunc(L[i] / std::pow(n, 2.0 / 21.0)) + 1;
      NumericVector x1 = seqNum(0, L1 - 1) / L1;
      NumericVector w1 = 1.0 - x1;
      NumericVector gai = seqNum(0, L1 - 1) * ga[Range(0, L1 - 1)] * w1;
      c2 = 3.0 * (2.0 * Rcpp::sum(Rcpp::pow(gai[Range(0, L1 - 1)], 2.0))) / (2.0 * M_PI);
      L[i + 1] = std::trunc(std::pow(n, 1.0 / 3.0) * std::pow(c2 / c1, 1.0 / 3.0)) + 1;

      if (L[i + 1] == L[i]) {
        runc = 0;
        LGopt = L[i + 1];
      }


    }
  }
  if (runc == 1) {
    LGopt = L[nit];
  }

  L1 = std::trunc(LGopt / std::pow(n, 2.0 / 21.0)) + 1;
  NumericVector x1 = seqNum(0, L1 - 1) / L1;
  NumericVector w1 = 1.0 - x1;
  NumericVector ga1 = seqNum(0, L1 - 1) * ga[Range(0, L1 - 1)] * w1;
  c20 = 3.0 * std::pow(2.0 * Rcpp::sum((ga1[Range(0, L1 - 1)])), 2.0) / (2.0 * M_PI);
  NumericVector w0 = (1.0 + cos(M_PI * x1)) / 2.0;
  NumericVector ga0 = ga[Range(0, L1 - 1)] * w0;
  c10 = std::pow(2.0 * Rcpp::sum((ga0[Range(0, L1 - 1)])) - ga0[0], 2.0) / (2.0 * M_PI);
  int L0opt = std::trunc(std::pow(n, 1.0 / 3.0) * std::pow((c20 / c10) / 2.0, 1.0 / 3.0)) + 1;
  NumericVector wacf = Rcpp::rev(seqNum(1, L0opt + 1)) / (L0opt + 1.0);
  NumericVector acfX = acovCpp(Xt, L0opt);
  double cf0LW = 2.0 * sum(acfX * wacf) - acfX[0];
  List out = List::create(_["cf0.LW"] = cf0LW, _["L0.opt"] = L0opt,
    _["LG.opt"] = LGopt);

  return out;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R

*/
