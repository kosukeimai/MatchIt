#include <Rcpp.h>
using namespace Rcpp;

// Functions for computing distance matrices. Slower than dist_to_matrix(dist(x))
// so don't bother using.

// [[Rcpp::export]]
NumericMatrix distC(NumericMatrix x){
  int p = x.ncol();
  int n = x.nrow();
  NumericMatrix out(n, n);
  NumericVector v1(p), v2(p);
  int i, j, k;
  double d, d0;

  for (i = 0; i < n; i++){
    v1 = x.row(i);
    for (j = i + 1; j < n; j++) {
      v2 = x.row(j);
      d = 0.0;
      for (k = 0; k < p; k++) {
        d0 = v1[k] - v2[k];
        d += d0 * d0;
      }
      d = sqrt(d);
      out(i, j) = d;
      out(j, i) = d;
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix distC2(NumericMatrix x, NumericMatrix y){
  int p = x.ncol();
  NumericMatrix out(x.nrow(), y.nrow());
  NumericVector v1(p), v2(p);
  int i, j, k;
  double d, d0;

  for (i = 0; i < x.nrow(); i++){
    v1 = x.row(i);
    for (j = 0; j < y.nrow(); j++) {
      v2 = y.row(j);
      d = 0.0;
      for (k = 0; k < p; k++) {
        d0 = v1[k] - v2[k];
        d += d0 * d0;
      }
      out(i, j) = sqrt(d);
    }
  }
  return out;
}