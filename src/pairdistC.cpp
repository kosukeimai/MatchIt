#include "internal.h"
using namespace Rcpp;

// [[Rcpp::export]]
double pairdistsubC(const NumericVector& x,
                    const IntegerVector& t,
                    const IntegerVector& s) {

  double dist = 0;

  R_xlen_t i, j;
  int s_i, o_i, o_j;
  int k = 0;

  //`base::order()`'s radix sort beats every C++ alternative measured here by 3-8x at
  //these sizes; see _dev/cpp-cleanup-notes.md. Looked up in the base environment
  //because `Function("order")` searches from the global environment, where a user
  //object of that name would mask it.
  Function ord = Environment::base_env()["order"];
  IntegerVector o = ord(s);
  o = o - 1;

  R_xlen_t n = sum(!is_na(s));

  for (i = 0; i < n; i++) {
    o_i = o[i];
    s_i = s[o_i];

    for (j = i + 1; j < n; j++) {
      o_j = o[j];

      if (s[o_j] != s_i) {
        break;
      }

      if (t[o_j] == t[o_i]) {
        continue;
      }

      //Numerically stable formula for adding new observation to a mean
      k++;
      dist += (std::abs(x[o_j] - x[o_i]) - dist) / k;
    }
  }

  return dist;
}
