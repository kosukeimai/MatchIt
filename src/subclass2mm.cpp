#include <Rcpp.h>
#include "internal.h"
using namespace Rcpp;

//Turns subclass vector given as a factor into a numeric match.matrix.
//focal is the treatment level (0/1) that corresponds to the rownames.

// [[Rcpp::export]]
IntegerMatrix subclass2mmC(const IntegerVector& subclass_,
                           const IntegerVector& treat,
                           const int& focal) {

  LogicalVector na_sub = is_na(subclass_);
  IntegerVector unique_sub = unique(as<IntegerVector>(subclass_[!na_sub]));
  IntegerVector subclass = match(subclass_, unique_sub) - 1;

  int nsub = unique_sub.size();

  int n = treat.size();
  IntegerVector ind = Range(0, n - 1);
  IntegerVector ind_focal = ind[treat == focal];
  int n1 = ind_focal.size();

  IntegerVector subtab = rep(-1, nsub);

  int i;
  for (i = 0; i < n; i++) {
    if (na_sub[i]) {
      continue;
    }

    subtab[subclass[i]]++;
  }

  int mm_col = max(subtab);

  IntegerMatrix mm(n1, mm_col);
  mm.fill(NA_INTEGER);
  CharacterVector lab = treat.names();

  IntegerVector ss = rep(NA_INTEGER, n1);

  int s, si;
  for (i = 0; i < n1; i++) {
    if (na_sub[ind_focal[i]]) {
      continue;
    }

      ss[i] = subclass[ind_focal[i]];
  }

  for (i = 0; i < n; i++) {
    if (treat[i] == focal) {
      continue;
    }

    if (na_sub[i]) {
      continue;
    }

    si = subclass[i];

    for (s = 0; s < n1; s++) {
      if (!std::isfinite(ss[s])) {
        continue;
      }

      if (si != ss[s]) {
        continue;
      }

      mm(s, sum(!is_na(mm(s, _)))) = i;
      break;
    }
  }

  mm = mm + 1;

  rownames(mm) = lab[ind_focal];

  return mm;
}