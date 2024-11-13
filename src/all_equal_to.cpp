#include "internal.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// Templated function to check if a vector has exactly n unique values
template <int RTYPE>
bool all_equal_to_(Vector<RTYPE> x,
                   typename traits::storage_type<RTYPE>::type y) {

  for (auto xi : x) {
    if (xi != y) {
      return false;
    }
  }

  return true;
}

// Wrapper function to handle different types of R vectors
// [[Rcpp::export]]
bool all_equal_to(RObject x,
                  RObject y) {

  switch (TYPEOF(x)) {
  case INTSXP:
    return all_equal_to_<INTSXP>(as<IntegerVector>(x), as<int>(y));
  case REALSXP:
    return all_equal_to_<REALSXP>(as<NumericVector>(x), as<double>(y));
  case LGLSXP:
    return all_equal_to_<LGLSXP>(as<LogicalVector>(x), as<bool>(y));
  default:
    stop("Unsupported vector type");
  }
}