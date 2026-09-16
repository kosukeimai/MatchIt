#include "internal.h"
using namespace Rcpp;

//Turns subclass vector given as a factor into a numeric match.matrix.
//focal is the treatment level (0/1) that corresponds to the rownames.

// [[Rcpp::export]]
IntegerMatrix subclass2mmC(const IntegerVector& subclass_,
                           const IntegerVector& treat,
                           int focal) {

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
  const CharacterVector lab = treat.names();

  //First row of `mm` belonging to each subclass. The rows are scanned in order and
  //the first match wins, which is what the loop this replaces did with its `break`;
  //doing it once makes the assignment below O(n) instead of O(n * n1).
  std::vector<int> first_row(nsub, -1);

  for (i = 0; i < n1; i++) {
    if (na_sub[ind_focal[i]]) {
      continue;
    }

    int si = subclass[ind_focal[i]];

    if (first_row[si] < 0) {
      first_row[si] = i;
    }
  }

  //Next column to fill in each row, rather than recomputing it with
  //`sum(!is_na(mm(s, _)))`, which allocates
  std::vector<int> mm_filled(n1, 0);

  for (i = 0; i < n; i++) {
    if (treat[i] == focal) {
      continue;
    }

    if (na_sub[i]) {
      continue;
    }

    int s = first_row[subclass[i]];

    if (s < 0) {
      continue;
    }

    mm(s, mm_filled[s]++) = i;
  }

  mm = mm + 1;
  rownames(mm) = as<CharacterVector>(lab[ind_focal]);

  return mm;
}

// [[Rcpp::export]]
IntegerVector mm2subclassC(const IntegerMatrix& mm,
                           const IntegerVector& treat,
                           const Nullable<int>& focal = R_NilValue) {

  const CharacterVector lab = treat.names();

  R_xlen_t n1 = treat.size();

  IntegerVector subclass(n1);
  subclass.fill(NA_INTEGER);
  subclass.names() = lab;

  const IntegerVector ind1 = focal.isNotNull() ?
  which(treat == as<int>(focal)) :
    IntegerVector(match(as<CharacterVector>(rownames(mm)), lab) - 1);

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
