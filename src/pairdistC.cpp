#include <Rcpp.h>
#include "internal.h"
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
  LogicalVector in_s_i(n);
  NumericVector x_t0(n);
  IntegerVector t_ind_s(n), c_ind_s(n);

  int k = 0;
  int i, i1, n1_s;
  for (i = 1; i <= num_sub; ++i) {
    in_s_i = (s == i);

    t_ind_s = which(t == 1 & in_s_i);
    c_ind_s = which(t == 0 & in_s_i);

    n1_s = t_ind_s.size();

    x_t0 = x[c_ind_s];

    for (i1 = 0; i1 < n1_s; ++i1) {
      dist += sum(Rcpp::abs(x[t_ind_s[i1]] - x_t0));
    }
    k += n1_s * c_ind_s.size();
  }

  dist /= k;

  return dist;
}