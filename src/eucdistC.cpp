#include "internal.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector eucdistC_N1xN0(const NumericMatrix& x,
                             const IntegerVector& t) {

  IntegerVector ind0 = which(t == 0);
  IntegerVector ind1 = which(t == 1);
  int p = x.ncol();

  NumericVector dist(ind1.size() * ind0.size());

  R_xlen_t k = 0;
  for (int i0 : ind0) {
    for (int i1 : ind1) {
      double d = 0;
      for (int i = 0; i < p; i++) {
        double di = x(i0, i) - x(i1, i);
        d += di * di;
      }
      dist[k] = std::sqrt(d);
      k++;
    }
  }

  dist.attr("dim") = Dimension(ind1.size(), ind0.size());

  return dist;
}
