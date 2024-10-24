#include <Rcpp.h>
#include "internal.h"
using namespace Rcpp;

// [[Rcpp::export]]
double pairdistsubC(const NumericVector& x,
                    const IntegerVector& t,
                    const IntegerVector& s) {

  double dist = 0;

  R_xlen_t i, j;
  int s_i, ord_i, ord_j;
  int k = 0;

  Function o("order");
  IntegerVector ord = o(s);
  ord = ord - 1;

  R_xlen_t n = sum(!is_na(s));


  for (i = 0; i < n; i++) {
    ord_i = ord[i];
    s_i = s[ord_i];

    for (j = i + 1; j < n; j++) {
      ord_j = ord[j];

      if (s[ord_j] != s_i) {
        break;
      }

      if (t[ord_j] == t[ord_i]) {
        continue;
      }

      k++;
      dist += (std::abs(x[ord_j] - x[ord_i]) - dist) / k;
    }
  }

  return dist;
}