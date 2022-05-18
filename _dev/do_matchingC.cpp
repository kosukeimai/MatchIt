#include <Rcpp.h>
#include "internal.h"
using namespace Rcpp;

//Performs nearest neighbor matching given an n1 by n0 distance matrix,
//with Inf values in the matrix forbidding matches.
//
//Arguments:
//
// distmat - n1xn0 numeric distance matrix
// treat - integer vector of 0/1 treatment status
// ord - integer vector of length n1 indicating matching order for treated units
// reuse_max - number of times each control unit can be matched set to >=n1 for
//             matching with replacement (unrestricted)
// ratio - integer vector of length n1 with number of controls for each treated
//         unit (is a vector to enable variable ratio matching)

// [[Rcpp::export]]
IntegerMatrix do_matchingC(const NumericMatrix& distmat,
                            const IntegerVector& treat,
                            const IntegerVector& ord,
                            const int& reuse_max,
                            const IntegerVector& ratio) {

  int max_rat = max(ratio);
  int n1 = sum(treat == 1);
  int n0 = sum(treat == 0);
  IntegerMatrix mm(n1, max_rat);
  mm.fill(NA_INTEGER);
  rownames(mm) = rownames(distmat);

  NumericMatrix distmat_ = clone(distmat);

  IntegerVector num_matches_t(n1);
  IntegerVector num_matches_c(n0);
  IntegerVector ind0_in_sample = which(treat == 0); //Index of control units

  int r, i, i_, m;
  double dm;

  for (r = 0; r < max_rat; r++) {
    for (i_ = 0; i_ < n1; i_++) {

      i = ord[i_] - 1; //-1 because it is in R indices

      //If treated units has received enough matched, move on
      if (num_matches_t[i] == ratio[i]) continue;

      //Find closest control index
      m = which_min(distmat_.row(i));

      dm = distmat_(i, m);

      //If closest is Inf, move on and prevent future attempts
      if (Rcpp::traits::is_infinite<REALSXP>(dm)) {
        num_matches_t[i] = ratio[i]; //prevents future matches
        continue;
      }

      //Assign matched index to matchmatrix
      mm(i, r) = ind0_in_sample[m] + 1; //+1 because we want to return R indices

      //Prevent re-matching by setting distance to Inf
      distmat_(i, m) = R_PosInf;

      //Increment number of times treated and matched control have been matched
      num_matches_t[i] = num_matches_t[i] + 1;
      num_matches_c[m] = num_matches_c[m] + 1;

      //If control units has been matched enough times, prevent future matches
      //by setting all distances to Inf
      if (num_matches_c[m] == reuse_max) {
        distmat_.column(m) = rep(R_PosInf, n1);
      }
    }
  }

  return mm;
}