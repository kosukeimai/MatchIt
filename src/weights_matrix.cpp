#include <Rcpp.h>
using namespace Rcpp;

// Computes matching weights from match.matrix

// [[Rcpp::export]]
NumericVector weights_matrix(const IntegerMatrix& mm,
                             const IntegerVector& treat) {
  int n = treat.size();
  IntegerVector ind = Range(0, n - 1);
  IntegerVector ind0 = ind[treat == 0];
  IntegerVector ind1 = ind[treat == 1];

  NumericVector weights (n);
  weights.fill(0);

  int nr = mm.nrow();
  int nc = mm.ncol();

  int r, c, row_not_na, which_c, t_ind;
  double weights_c, add_w;
  IntegerVector row_r(nc);

  for (r = 0; r < nr; r++) {
    row_r = na_omit(mm(r, _));
    row_not_na = row_r.size();
    add_w = 1.0/static_cast<double>(row_not_na);
    if (row_not_na == 0) {
      continue;
    }

    for (c = 0; c < row_not_na; c++) {
      which_c = row_r[c] - 1;
      weights_c = weights[which_c];
      weights[which_c] = weights_c + add_w;
    }

    t_ind = ind1[r];
    weights[t_ind] = 1;
  }

  NumericVector c_weights = weights[ind0];
  double sum_c_w = sum(c_weights);
  double sum_matched_c = sum(c_weights > 0);
  int n0 = ind0.size();

  for (int i = 0; i < n0; i++ ) {
    which_c = ind0[i];
    weights[which_c] = c_weights[i] * sum_matched_c / sum_c_w;
  }

  return weights;
}