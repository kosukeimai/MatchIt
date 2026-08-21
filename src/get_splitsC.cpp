#include "internal.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector get_splitsC(const NumericVector& x,
                          double caliper) {

  NumericVector x_ = unique(x);
  NumericVector x_sorted = x_.sort();

  R_xlen_t n = x_sorted.size();

  if (n <= 1) {
    return NumericVector(0);
  }

  //Accumulated in a std::vector because NumericVector::push_back() reallocates and
  //copies the whole vector on every call, which is quadratic in the number of splits
  std::vector<double> splits;

  splits.push_back(x_sorted[0]);

  for (R_xlen_t i = 1; i < n; i++) {
    if (x_sorted[i] - x_sorted[i - 1] <= caliper) continue;

    splits.push_back((x_sorted[i] + x_sorted[i - 1]) / 2);
  }

  splits.push_back(x_sorted[n - 1]);

  return wrap(splits);
}
