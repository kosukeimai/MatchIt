#include <Rcpp.h>
#include "internal.h"
#include <cmath>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// Computes matching weights from subclass
// [[Rcpp::export]]
NumericVector weights_subclassC(const IntegerVector& subclass_,
                                const IntegerVector& treat_,
                                const Nullable<int>& focal_ = R_NilValue) {

  CharacterVector lab = treat_.names();
  IntegerVector unique_treat = {0, 1};
  // IntegerVector unique_treat = unique(treat_);
  // std::sort(unique_treat.begin(), unique_treat.end());
  int g = unique_treat.size();
  // IntegerVector treat = match(treat_, unique_treat) - 1;
  IntegerVector treat = treat_;

  IntegerVector non_na_sub = which(!is_na(subclass_));
  // IntegerVector unique_sub = unique(as<IntegerVector>(subclass_[non_na_sub]));
  // IntegerVector subclass = match(subclass_, unique_sub) - 1;
  IntegerVector subclass = subclass_ - 1;

  int nsub = max(as<IntegerVector>(subclass_[non_na_sub]));

  R_xlen_t n = treat.size();
  int gi;

  NumericVector weights(n);
  weights.fill(0.0);
  weights.names() = lab;

  double subtab[nsub][g];
  NumericMatrix subtab(nsub, g);
  subtab.fill(0.0);

  //Count number of units in each subclass by treatment
  for (int i : non_na_sub) {
    subtab(subclass[i], treat[i])++;
  }

  NumericVector subtab_total = rowSums(subtab);

  NumericMatrix subtab(nsub, g);
  subtab.fill(0.0);

  NumericVector weights_gi;
  IntegerVector indg;
  double sum_w;
  double sum_matched;

  if (focal_.isNull()) {
    for (gi = 0; gi < g; gi++) {
      sub_weight(_, gi) = subtab_total / subtab(_, gi);
    }

    for (int i : non_na_sub) {
      weights[i] = sub_weight(subclass[i], treat[i]);
    }
  }
  else {
    int focal = as<int>(focal_);
    sub_weight.fill(1);

    for (gi = 0; gi < g; gi++) {
      if (gi != focal) {
        sub_weight(_, gi) = subtab(_, focal) / subtab(_, gi);
      }
    }

    for (int i : non_na_sub) {
      weights[i] = sub_weight(subclass[i], treat[i]);
    }
  }

  return weights;
}