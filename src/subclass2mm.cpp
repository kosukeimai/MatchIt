#include "internal.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

//Turns subclass vector given as a factor into a numeric match.matrix.
//focal is the treatment level (0/1) that corresponds to the rownames.

// [[Rcpp::export]]
IntegerMatrix subclass2mmC(const IntegerVector& subclass_,
                           const IntegerVector& treat,
                           const int& focal) {

  LogicalVector na_sub = is_na(subclass_);
  IntegerVector unique_sub = unique(as<IntegerVector>(subclass_[!na_sub]));
  IntegerVector subclass = match(subclass_, unique_sub) - 1;

  R_xlen_t nsub = unique_sub.size();

  R_xlen_t n = treat.size();
  IntegerVector ind = Range(0, n - 1);
  IntegerVector ind_focal = ind[treat == focal];
  R_xlen_t n1 = ind_focal.size();

  IntegerVector subtab(nsub);
  subtab.fill(-1);

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
  rownames(mm) = as<CharacterVector>(lab[ind_focal]);

  return mm;
}

// [[Rcpp::export]]
IntegerVector mm2subclassC(const IntegerMatrix& mm,
                           const IntegerVector& treat,
                           const Nullable<int>& focal = R_NilValue) {

  CharacterVector lab = treat.names();

  R_xlen_t n1 = treat.size();

  IntegerVector subclass(n1);
  subclass.fill(NA_INTEGER);
  subclass.names() = lab;

  IntegerVector ind1;
  if (focal.isNotNull()) {
    ind1 = which(treat == as<int>(focal));
  }
  else {
    ind1 = match(as<CharacterVector>(rownames(mm)), lab) - 1;
  }

  R_xlen_t r = mm.nrow();
  R_xlen_t ki = 0;
  int ri;

  IntegerVector s(r);
  std::vector<std::string> levs;
  levs.reserve(r);

  for (R_xlen_t i : which(!is_na(mm))) {
    ri = i % r; //row

    //If first in column, assign subclass
    if (i / r == 0) {
      ki++;

      s[ri] = ki;
      subclass[ind1[ri]] = ki;

      levs.push_back(std::to_string(ki));
    }

    subclass[mm[i] - 1] = s[ri];
  }

  subclass.attr("class") = "factor";
  subclass.attr("levels") = levs;

  return subclass;
}