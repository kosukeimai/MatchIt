#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// Rcpp internal functions

//C implementation of tabulate(). Faster than base::tabulate(), but real
//use is in subclass2mmC().

// [[Rcpp::interfaces(cpp)]]
IntegerVector tabulateC_(const IntegerVector& bins,
                         const Nullable<int>& nbins = R_NilValue) {
  int max_bin;

  if (nbins.isNotNull()) max_bin = as<int>(nbins);
  else max_bin = max(na_omit(bins));

  IntegerVector counts(max_bin);
  int n = bins.size();
  for (int i = 0; i < n; i++) {
    if (bins[i] > 0 && bins[i] <= max_bin) {
      counts[bins[i] - 1]++;
    }
  }

  return counts;
}

//Rcpp port of base::which

// [[Rcpp::interfaces(cpp)]]
IntegerVector which(const LogicalVector& x) {
  IntegerVector ind = Range(0, x.size() - 1);
  return ind[x];
}

// [[Rcpp::interfaces(cpp)]]
bool antiexact_okay(const int& aenc,
                    const int& ii,
                    const int& i,
                    const IntegerMatrix& antiexact_covs) {
  if (aenc == 0) {
    return true;
  }

  for (int j = 0; j < aenc; j++) {
    if (antiexact_covs(ii, j) == antiexact_covs(i, j)) {
      return false;
    }
  }

  return true;
}

// [[Rcpp::interfaces(cpp)]]
bool caliper_covs_okay(const int& ncc,
                       const int& ii,
                       const int& i,
                       const NumericMatrix& caliper_covs_mat,
                       const NumericVector& caliper_covs) {
  if (ncc == 0) {
    return true;
  }

  for (int j = 0; j < ncc; j++) {
    if (caliper_covs[j] >= 0) {
      if (std::abs(caliper_covs_mat(ii, j) - caliper_covs_mat(i, j)) > caliper_covs[j]) {
        return false;
      }
    }
    else {
      if (std::abs(caliper_covs_mat(ii, j) - caliper_covs_mat(i, j)) <= -caliper_covs[j]) {
        return false;
      }
    }
  }

  return true;
}

// [[Rcpp::interfaces(cpp)]]
bool mm_okay(const int& r,
             const int& i,
             const IntegerVector& mm_ordi) {

  if (r > 1) {
    for (int j : mm_ordi) {
      if (i == j) {
        return false;
      }
    }
  }

  return true;
}

// [[Rcpp::interfaces(cpp)]]
bool exact_okay(const bool& use_exact,
                const int& ii,
                const int& i,
                const IntegerVector& exact) {

  if (!use_exact) {
    return true;
  }

  return exact[ii] == exact[i];
}

// [[Rcpp::interfaces(cpp)]]
int find_both(const int& t_id,
              const IntegerVector& ind_d_ord,
              const IntegerVector& match_d_ord,
              const IntegerVector& treat,
              const NumericVector& distance,
              const LogicalVector& eligible,
              const int& gi,
              const int& r,
              const IntegerVector& mm_ordi,
              const int& ncc,
              const NumericMatrix& caliper_covs_mat,
              const NumericVector& caliper_covs,
              const double& caliper_dist,
              const bool& use_exact,
              const IntegerVector& exact,
              const int& aenc,
              const IntegerMatrix& antiexact_covs,
              const int& first_control,
              const int& last_control) {

  int ii = match_d_ord[t_id];

  int min_ii = first_control;
  int max_ii = last_control;

  int iil = ii;
  int iir = ii;
  int il = -1;
  int ir = -1;

  bool l_found = (iil <= min_ii);
  bool r_found = (iir >= max_ii);

  double di = distance[t_id];
  double distl, distr;

  while (!l_found || !r_found) {
    if (!l_found) {
      if (iil == min_ii) {
        l_found = true;
        il = -1;
      }
      else {
        iil--;
        il = ind_d_ord[iil];

        //Left
        if (eligible[il]) {
          if (treat[il] == gi) {
            if (mm_okay(r, il, mm_ordi)) {

              distl = std::abs(di - distance[il]);

              if (r_found && ir >= 0 && distl > distr) {
                return ir;
              }

              if (distl > caliper_dist) {
                il = -1;
                l_found = true;
              }
              else {
                if (exact_okay(use_exact, t_id, il, exact)) {
                  if (antiexact_okay(aenc, t_id, il, antiexact_covs)) {
                    if (caliper_covs_okay(ncc, t_id, il, caliper_covs_mat, caliper_covs)) {
                      l_found = true;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    if (!r_found) {
      if (iir == max_ii) {
        r_found = true;
        ir = -1;
      }
      else {
        iir++;
        ir = ind_d_ord[iir];

        //Right
        if (eligible[ir]) {
          if (treat[ir] == gi) {
            if (mm_okay(r, ir, mm_ordi)) {

              distr = std::abs(di - distance[ir]);

              if (l_found && il >= 0 && distl <= distr) {
                return il;
              }

              if (distr > caliper_dist) {
                ir = -1;
                r_found = true;
              }
              else {
                if (exact_okay(use_exact, t_id, ir, exact)) {
                  if (antiexact_okay(aenc, t_id, ir, antiexact_covs)) {
                    if (caliper_covs_okay(ncc, t_id, ir, caliper_covs_mat, caliper_covs)) {
                      r_found = true;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  if (il < 0) {
    return ir;
  }

  if (ir < 0) {
    return il;
  }

  if (distl <= distr) {
    return il;
  }
  else {
    return ir;
  }
}

// [[Rcpp::interfaces(cpp)]]
int find_lr(const int& prev_match,
            const int& t_id,
            const IntegerVector& ind_d_ord,
            const IntegerVector& match_d_ord,
            const IntegerVector& treat,
            const NumericVector& distance,
            const LogicalVector& eligible,
            const int& gi,
            const int& ncc,
            const NumericMatrix& caliper_covs_mat,
            const NumericVector& caliper_covs,
            const double& caliper_dist,
            const bool& use_exact,
            const IntegerVector& exact,
            const int& aenc,
            const IntegerMatrix& antiexact_covs,
            const int& first_control,
            const int& last_control) {

  int ik, iik;
  double dist;

  int prev_pos;
  int ii = match_d_ord[t_id];

  int z;
  if (prev_match < 0) {
    if (prev_match == -1) {
      z = -1;
    }
    else {
      z = 1;
    }

    prev_pos = ii + z;
  }
  else {
    prev_pos = match_d_ord[prev_match];
    if (prev_pos < ii) {
      z = -1;
    }
    else {
      z = 1;
    }
  }

  int min_ii = 0;
  int max_ii = ind_d_ord.size() - 1;

  if (z == -1) {
    min_ii = first_control;
  }
  else {
    max_ii = last_control;
  }

  for (iik = prev_pos; iik >= min_ii && iik <= max_ii; iik = iik + z) {
    ik = ind_d_ord[iik];

    if (!eligible[ik]) {
      continue;
    }

    if (treat[ik] != gi) {
      continue;
    }

    dist = std::abs(distance[t_id] - distance[ik]);

    if (dist > caliper_dist) {
      return -1;
    }

    if (!exact_okay(use_exact, t_id, ik, exact)) {
      continue;
    }

    if (!antiexact_okay(aenc, t_id, ik, antiexact_covs)) {
      continue;
    }

    if (!caliper_covs_okay(ncc, t_id, ik, caliper_covs_mat, caliper_covs)) {
      continue;
    }

    return ik;
  }

  return -1;
}

// [[Rcpp::interfaces(cpp)]]
IntegerVector swap_pos(IntegerVector x,
                       const int& a,
                       const int& b) {
  int xa = x[a];

  x[a] = x[b];
  x[b] = xa;

  return x;
}

// [[Rcpp::interfaces(cpp)]]
double max_finite(const NumericVector& x) {
  double m = NA_REAL;

  int n = x.size();
  int i;
  bool found = false;

  //Find first finite value
  for (i = 0; i < n; i++) {
    if (std::isfinite(x[i])) {
      m = x[i];
      found = true;
      break;
    }
  }

  //If none found, return NA
  if (!found) {
    return m;
  }

  //Find largest finite value
  for (int j = i + 1; j < n; j++) {
    if (!std::isfinite(x[i])) {
      continue;
    }

    if (x[j] > m) {
      m = x[j];
    }
  }

  return m;
}

// [[Rcpp::interfaces(cpp)]]
double min_finite(const NumericVector& x) {
  double m = NA_REAL;

  int n = x.size();
  int i;
  bool found = false;

  //Find first finite value
  for (i = 0; i < n; i++) {
    if (std::isfinite(x[i])) {
      m = x[i];
      found = true;
      break;
    }
  }

  //If none found, return NA
  if (!found) {
    return m;
  }

  //Find smallest finite value
  for (int j = i + 1; j < n; j++) {
    if (!std::isfinite(x[i])) {
      continue;
    }

    if (x[j] < m) {
      m = x[j];
    }
  }

  return m;
}