#include "internal.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// Templated function to check if a vector has exactly n unique values
template <int RTYPE>
bool has_n_unique_(Vector<RTYPE> x,
                   const int& n) {

  Vector<RTYPE> seen(n);
  seen[0] = x[0];
  int n_seen = 1;

  int j;
  bool was_seen;

  // Iterate over the vector and add elements to the unordered set
  for (auto it = x.begin() + 1; it != x.end(); ++it) {
    if (*it == *(it - 1)) {
      continue;
    }

    was_seen = false;

    for (j = 0; j < n_seen; j++) {
      if (*it == seen[j]) {
        was_seen = true;
        break;
      }
    }

    if (!was_seen) {
      n_seen++;

      if (n_seen > n) {
        return false;
      }

      seen[n_seen - 1] = *it;
    }
  }

  // Check if the number of unique elements is exactly n
  return n_seen == n;
}

// Wrapper function to handle different types of R vectors
// [[Rcpp::export]]
bool has_n_unique(const SEXP& x,
                  const int& n) {
  switch (TYPEOF(x)) {
  case INTSXP:
    return has_n_unique_<INTSXP>(x, n);
  case REALSXP:
    return has_n_unique_<REALSXP>(x, n);
  case STRSXP:
    return has_n_unique_<STRSXP>(x, n);
  case LGLSXP:
    return has_n_unique_<LGLSXP>(x, n);
  default:
    stop("Unsupported vector type");
  }
}