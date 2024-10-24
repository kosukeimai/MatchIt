#include <Rcpp.h>
#include "internal.h"
#include <cmath>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// Computes matching weights from match.matrix
// [[Rcpp::export]]
NumericVector weights_matrixC(const IntegerMatrix& mm,
                              const IntegerVector& treat_,
                              const Nullable<int>& focal = R_NilValue) {

  CharacterVector lab = treat_.names();
  IntegerVector unique_treat = unique(treat_);
  std::sort(unique_treat.begin(), unique_treat.end());
  int g = unique_treat.size();
  IntegerVector treat = match(treat_, unique_treat) - 1;

  R_xlen_t n = treat.size();
  int gi;

  NumericVector weights(n);
  weights.fill(0.0);
  weights.names() = lab;

  IntegerVector row_ind;
  if (focal.isNotNull()) {
    row_ind = which(treat == as<int>(focal));
  }
  else {
    row_ind = match(as<CharacterVector>(rownames(mm)), lab) - 1;
  }

  NumericVector matches_g = rep(0.0, g);

  IntegerVector row_r(mm.ncol());

  for (int r : which(!is_na(mm(_, 0)))) {

    row_r = na_omit(mm.row(r));

    for (gi = 0; gi < g; gi++) {
      matches_g[gi] = 0.0;
    }

    for (int i : row_r - 1) {
      matches_g[treat[i]] += 1.0;
    }

    for (int i : row_r - 1) {
      if (matches_g[treat[i]] == 0.0) {
        continue;
      }

      weights[i] += 1.0/matches_g[treat[i]];
    }

    weights[row_ind[r]] += 1.0;
  }

  //Scale control weights to sum to number of matched controls
  NumericVector weights_gi;
  IntegerVector indg;
  double sum_w;
  double sum_matched;

  for (gi = 0; gi < g; gi++) {
    indg = which(treat == gi);

    weights_gi = weights[indg];

    sum_w = sum(weights_gi);

    if (sum_w == 0) {
      continue;
    }

    sum_matched = sum(weights_gi > 0);

    if (sum_matched == sum_w) {
      continue;
    }

    for (int i : indg) {
      weights[i] *= sum_matched / sum_w;
    }
  }

  return weights;
}