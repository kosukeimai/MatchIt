#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double pairdistsubC(const NumericVector& x_,
                    const IntegerVector& t_,
                    const IntegerVector& s_,
                    const int& num_sub) {

  double dist = 0;

  LogicalVector not_na_sub = !is_na(s_);
  NumericVector x = x_[not_na_sub];
  IntegerVector t = t_[not_na_sub];
  IntegerVector s = s_[not_na_sub];

  int n = t.size();
  LogicalVector t_i(n), in_s_i(n);
  NumericVector x_t1(n), x_t0(n);

  int k = 0;
  int i, i1, i0;
  for (i = 1; i <= num_sub; ++i) {
    in_s_i = (s == i);

    t_i = (t == 1);
    t_i[!in_s_i] = false;
    x_t1 = x[t_i];

    t_i = (t == 0);
    t_i[!in_s_i] = false;
    x_t0 = x[t_i];

    for (i1 = 0; i1 < x_t1.size(); ++i1) {
      for (i0 = 0; i0 < x_t0.size(); ++i0) {
        dist += abs(x_t1[i1] - x_t0[i0]);
        ++k;
      }
    }
  }

  dist /= k;

  return dist;
}