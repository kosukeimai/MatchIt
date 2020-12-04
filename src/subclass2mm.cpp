#include <Rcpp.h>
#include <tabulateC.cpp>
using namespace Rcpp;

//Turns subclass vector given as a factor into a numeric match.matrix.
//focal is the treatment level (0/1) that corresponds to the rownames.

// [[Rcpp::export]]
IntegerMatrix subclass2mmC(const IntegerVector& subclass,
                           const IntegerVector& treat,
                           const int& focal) {

  IntegerVector tab = tabulateC(subclass);
  int mm_col = max(tab) - 1;

  IntegerVector ind = Range(0, treat.size() - 1);
  IntegerVector ind1 = ind[treat == focal];
  int n1 = ind1.size();

  IntegerMatrix mm(n1, mm_col);
  mm.fill(NA_INTEGER);
  rownames(mm) = as<CharacterVector>(treat.names())[ind1];

  IntegerVector in_sub = ind[!is_na(subclass)];
  IntegerVector ind_in_sub = ind[in_sub];
  IntegerVector ind0_in_sub = ind_in_sub[as<IntegerVector>(treat[in_sub]) == 0];
  IntegerVector sub0_in_sub = subclass[ind0_in_sub];

  int i, t, s, nmc, mci;
  IntegerVector mc(mm_col);

  for (i = 0; i < n1; i++) {
    t = ind1[i];
    s = subclass[t];

    if (s != NA_INTEGER) {
      mc = ind0_in_sub[sub0_in_sub == s];
      nmc = mc.size();
      for (mci = 0; mci < nmc; mci++) {
        mm(i, mci) = mc[mci] + 1;
      }
    }
  }

  return mm;
}

