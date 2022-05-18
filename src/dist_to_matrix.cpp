#include <Rcpp.h>
using namespace Rcpp;
//Faster alternative to stats:::as.matrix.dist().

// [[Rcpp::export]]
NumericMatrix dist_to_matrixC(const NumericVector& d) {
  int n = d.attr("Size");

  NumericMatrix m(n,n);

  double dk;
  int i, j;
  int k = 0;

  for (i = 0; i < n; i++) {
    for (j = i + 1; j < n; j++) {
      dk = d[k];
      m(i, j) = dk;
      m(j, i) = dk;
      k++;
    }
  }

  if (d.hasAttribute("Labels")) {
    CharacterVector lab = d.attr("Labels");
    rownames(m) = lab;
    colnames(m) = lab;
  }
  return m;
}
