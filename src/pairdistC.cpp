#include <Rcpp.h>
#include "internal.h"
using namespace Rcpp;

// [[Rcpp::export]]
double pairdistsubC(const NumericVector& x,
                    const IntegerVector& t,
                    const IntegerVector& s) {

  double dist = 0;

  int n = t.size();
  int i, j;
  int k = 0;

  for (i = 0; i < n - 1; i++) {
    if (!std::isfinite(s[i])) {
      continue;
    }

    for (j = i + 1; j < n; j++) {
      if (s[i] != s[j]) {
        continue;
      }

      if (t[i] == t[j]) {
        continue;
      }

      dist += std::abs(x[i] - x[j]);
      k++;
    }
  }

  dist /= k;

  return dist;
}