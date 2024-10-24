#include <Rcpp.h>
#include "internal.h"
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector tabulateC(const IntegerVector& bins,
                        const Nullable<int>& nbins = R_NilValue) {

  int nbins_ = 0;
  if (nbins.isNotNull()) nbins_ = as<int>(nbins);

  return tabulateC_(bins, nbins_);
}