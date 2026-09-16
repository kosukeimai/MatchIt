#include "internal.h"
using namespace Rcpp;

// Computes matching weights from match.matrix
// [[Rcpp::export]]
NumericVector weights_matrixC(const IntegerMatrix& mm,
                              const IntegerVector& treat_,
                              const Nullable<int>& focal = R_NilValue) {

  const CharacterVector lab = treat_.names();
  IntegerVector unique_treat = unique(treat_);
  std::sort(unique_treat.begin(), unique_treat.end());
  int g = unique_treat.size();
  IntegerVector treat = match(treat_, unique_treat) - 1;

  R_xlen_t n = treat.size();
  int gi;

  NumericVector weights(n);
  weights.fill(0.0);
  weights.names() = lab;

  //`treat` has been recoded to 0..g-1, so `focal` must be recoded the same way
  const IntegerVector row_ind = focal.isNotNull() ?
  which(treat == recode_focal(as<int>(focal), unique_treat)) :
    IntegerVector(match(as<CharacterVector>(rownames(mm)), lab) - 1);

  std::vector<double> matches_g(g, 0.0);

  R_xlen_t mm_ncol = mm.ncol();

  std::vector<int> row_r;
  row_r.reserve(mm_ncol);

  for (int r : which(!is_na(mm(_, 0)))) {

    //Filled by hand rather than with `na_omit(mm.row(r))`, which allocates two
    //vectors for every row of the match matrix
    row_r.clear();

    for (R_xlen_t j = 0; j < mm_ncol; j++) {
      if (mm(r, j) != NA_INTEGER) {
        row_r.push_back(mm(r, j) - 1);
      }
    }

    for (gi = 0; gi < g; gi++) {
      matches_g[gi] = 0.0;
    }

    for (int i : row_r) {
      matches_g[treat[i]] += 1.0;
    }

    for (int i : row_r) {
      if (matches_g[treat[i]] == 0.0) {
        continue;
      }

      weights[i] += 1.0/matches_g[treat[i]];
    }

    weights[row_ind[r]] += 1.0;
  }

  return weights;
}
