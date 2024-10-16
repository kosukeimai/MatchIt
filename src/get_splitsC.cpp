#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector get_splitsC(const NumericVector& x,
                          const double& caliper) {

  NumericVector splits;

  NumericVector x_ = unique(x);
  NumericVector x_sorted = x_.sort();

  int n = x_sorted.size();

  if (n <= 1) {
    return splits;
  }

  splits = x_sorted[0];

  for (int i = 1; i < x_sorted.length(); i++) {
    if (x_sorted[i] - x_sorted[i - 1] <= caliper) continue;

    splits.push_back((x_sorted[i] + x_sorted[i - 1]) / 2);
  }

  splits.push_back(x_sorted[n - 1]);

  return splits;
}