#include <Rcpp.h>
using namespace Rcpp;

//C implementation of tabulate. Faster than base::tabulate(), but real
//use is in subclass2mmC().

// [[Rcpp::export]]
IntegerVector tabulateC(const IntegerVector& bins,
                        const Nullable<int>& nbins = R_NilValue) {
  int max_bin;
  if (nbins.isNotNull()) max_bin = as<int>(nbins);
  else max_bin = max(na_omit(bins));

  IntegerVector counts(max_bin);
  int n = bins.size();
  for (int i = 0; i < n; i++) {
    if (bins[i] > 0 && bins[i] <= max_bin)
      counts[bins[i] - 1]++;
  }
  return counts;
}