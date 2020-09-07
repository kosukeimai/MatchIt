#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double pairdistsubC(const NumericVector& x_,
                    const IntegerVector& t_,
                    const IntegerVector& s_,
                    const int& num_sub) {

  double dist = 0;

  LogicalVector na_sub = is_na(s_);
  NumericVector x = x_[!na_sub];
  IntegerVector t = t_[!na_sub];
  IntegerVector s = s_[!na_sub];

  int n = t.size();
  IntegerVector s_unique = Range(0, num_sub - 1);
  LogicalVector t_i(n), in_s_i(n);
  NumericVector x_t1(n), x_t0(n);

  int k = 0;
  for (int i = 1; i <= num_sub; ++i) {
    in_s_i = s == i;
    in_s_i[is_na(s)] = false;

    t_i = (t == 1);
    t_i[!in_s_i] = false;
    x_t1 = x[t_i];

    t_i = (t == 0);
    t_i[!in_s_i] = false;
    x_t0 = x[t_i];

    for (int i1 = 0; i1 < x_t1.size(); ++i1) {
      for (int i0 = 0; i0 < x_t0.size(); ++i0) {
        dist += abs(x_t1[i1] - x_t0[i0]);
        ++k;
      }
    }
  }

  dist /= k;

  return dist;
}