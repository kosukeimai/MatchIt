#include "internal.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
IntegerVector subclass_scootC(const IntegerVector& subclass_,
                              const IntegerVector& treat_,
                              const NumericVector& x_,
                              const int& min_n) {

  if (min_n == 0) {
    return subclass_;
  }

  int m, i, s, s2;
  int best_i;
  double best_x, score;
  R_xlen_t nt;

  LogicalVector na_sub = is_na(subclass_);

  IntegerVector subclass = subclass_[!na_sub];
  IntegerVector treat = treat_[!na_sub];
  NumericVector x = x_[!na_sub];

  R_xlen_t n = subclass.size();

  IntegerVector unique_sub = unique(subclass);
  std::sort(unique_sub.begin(), unique_sub.end());

  subclass = match(subclass, unique_sub) - 1;

  R_xlen_t nsub = unique_sub.size();

  NumericVector subtab(nsub);
  IntegerVector indt;
  bool left = false;

  IntegerVector ut = unique(treat);

  for (int t : ut) {
    indt = which(treat == t);
    nt = indt.size();

    //Tabulate
    subtab.fill(0.0);

    for (int i : indt) {
      subtab[subclass[i]]++;
    }

    for (m = 0; m < min_n; m++) {
      while (min(subtab) <= 0) {
        for (s = 0; s < nsub; s++) {
          if (subtab[s] == 0) {
            break;
          }
        }

        //Find which way to look for new member
        if (s == nsub - 1) {
          left = true;
        }
        else if (s == 0) {
          left = false;
        }
        else {
          score = 0.;

          for (s2 = 0; s2 < nsub; s2++) {
            if (subtab[s2] <= 1) {
              continue;
            }

            if (s2 == s) {
              continue;
            }

            score += (subtab[s2] - 1) / static_cast<double>(s2 - s);
          }

          left = (score <= 0);
        }

        //Find which subclass to take from (s2)
        if (left) {
          for (s2 = s - 1; s2 >= 0; s2--) {
            if (subtab[s2] > 0) {
              break;
            }
          }
        }
        else {
          for (s2 = s + 1; s2 < nsub; s2++) {
            if (subtab[s2] > 0) {
              break;
            }
          }
        }

        //Find unit with closest x in that subclass to take
        for (i = 0; i < nt; i++) {
          if (subclass[indt[i]] == s2) {
            best_i = i;
            best_x = x[indt[i]];
            break;
          }
        }

        for (i = best_i + 1; i < nt; i++) {
          if (subclass[indt[i]] != s2) {
            continue;
          }

          if (left) {
            if (x[indt[i]] < best_x) {
              continue;
            }
          }
          else {
            if (x[indt[i]] > best_x) {
              continue;
            }
          }

          best_i = i;
          best_x = x[indt[i]];
        }

        subclass[indt[best_i]] = s;
        subtab[s]++;
        subtab[s2]--;
      }

      for (s = 0; s < nsub; s++) {
        subtab[s]--;
      }
    }
  }

  for (i = 0; i < n; i++) {
    subclass[i] = unique_sub[subclass[i]];
  }

  IntegerVector sub_out(subclass_.size());
  sub_out.fill(NA_INTEGER);

  sub_out[!na_sub] = subclass;

  return sub_out;
}