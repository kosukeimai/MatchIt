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

  R_xlen_t n = treat.size();
  IntegerVector ind = Range(0, n - 1);
  IntegerVector ind_focal = ind[treat == focal];
  R_xlen_t n1 = ind_focal.size();

  IntegerVector subtab = rep(-1, nsub);

  R_xlen_t i;
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

  IntegerVector ss(n1);
  ss.fill(NA_INTEGER);

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

// [[Rcpp::export]]
IntegerVector mm2subclassC(const IntegerMatrix& mm,
                           const IntegerVector& treat,
                           const Nullable<int>& focal = R_NilValue) {

  CharacterVector lab = treat.names();

  IntegerVector subclass(treat.size());
  subclass.fill(NA_INTEGER);
  subclass.names() = lab;

  IntegerVector ind1;
  if (focal.isNotNull()) {
    ind1 = which(treat == as<int>(focal));
  }
  else {
    ind1 = match(as<CharacterVector>(rownames(mm)), lab) - 1;
  }

  int r = mm.nrow();
  int ki = 0;

  for (int i : which(!is_na(mm))) {
    if (i / r == 0) {
      //If first entry in row, increment ki and assign subclass of treated
      ki++;
      subclass[ind1[i % r]] = ki;
    }

    subclass[mm[i] - 1] = ki;
  }

  CharacterVector levs(ki);
  for (int j = 0; j < ki; j++){
    levs[j] = std::to_string(j + 1);
  }

  subclass.attr("class") = "factor";
  subclass.attr("levels") = levs;

  return subclass;
}