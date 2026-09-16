#include <Rcpp.h>
#include "internal.h"
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector tabulateC(const IntegerVector& bins,
                        const Nullable<int>& nbins = R_NilValue) {

  const int nbins_ = nbins.isNotNull() ? as<int>(nbins) : 0;

  return tabulateC_(bins, nbins_);
}
