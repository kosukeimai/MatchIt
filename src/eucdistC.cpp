#include <Rcpp.h>
#include "internal.h"
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector eucdistC_N1xN0(const NumericMatrix& x,
                             const IntegerVector& t) {

  IntegerVector ind0 = which(t == 0);
  IntegerVector ind1 = which(t == 1);
  int p = x.ncol();
  int i;
  double d, di;

  NumericVector dist(ind1.size() * ind0.size());

  int k = 0;
  for (double i0 : ind0) {
    for (double i1 : ind1) {
      d = 0;
      for (i = 0; i < p; i++) {
        di = x(i0, i) - x(i1, i);
        d += di * di;
      }
      dist[k] = sqrt(d);
      k++;
    }
  }

  dist.attr("dim") = Dimension(ind1.size(), ind0.size());

  return dist;
}