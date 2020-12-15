#include <Rcpp.h>
#include "internal.h"
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector tabulateC(const IntegerVector& bins,
                        const Nullable<int>& nbins = R_NilValue) {
  return tabulateC_(bins, nbins);
}