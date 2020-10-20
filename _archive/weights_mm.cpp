#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector weights_mm(CharacterMatrix mm,
                         IntegerVector treat) {
  int n = treat.size();
  int n1 = mm.nrow();
  CharacterVector mrownames = rownames(mm);
  NumericVector weights(n);
  weights.names() = treat.names();
  IntegerVector ind = Range(0, n-1);
  ind.names() = weights.names();

  CharacterVector row_i(mm.ncol());
  NumericVector nmatches_i(1);
  String mrowname_i;
  String matched_controls_j;
  CharacterVector matched_controls_i(mm.ncol());
  int i, j, k;

  for (i = 0; i < n1; i++) {
    row_i = mm(i , _);
    mrowname_i = mrownames[i];
    nmatches_i[1] = sum(!is_na(row_i));

    if (nmatches_i[1] == 0) continue;
    weights[mrowname_i] = 1;
    matched_controls_i = na_omit(row_i);

    for (j = 0; j < nmatches_i[1]; j++) {
      matched_controls_j = matched_controls_i[j];
      k = weights.findName(matched_controls_j);
      weights[k] += 1/nmatches_i[1];
    }
  }

  return weights;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
system.time({
  w <- weights_mm(m$match.matrix, m$treat)
})
system.time({
  w2 <- weights.matrix(m$match.matrix, m$treat)
})

*/
