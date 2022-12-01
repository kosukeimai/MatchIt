// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <Rcpp.h>
#include "internal.h"
using namespace Rcpp;

// Version of nn_matchC that works when `distance` is a vector.
// Doesn't accept `distance_mat_` or `mah_covs_`.

bool check_in(int x,
              IntegerVector table) {
  for (int j = 0; j < table.size(); j++) {
    if (x == table[j]) return true;
  }
  return false;
}

int find_right(int ii,
               int last_control,
               IntegerVector treat,
               LogicalVector can_be_matched,
               int r,
               IntegerVector row,
               IntegerVector d_ord,
               NumericVector distance,
               bool use_caliper_dist,
               double caliper_dist,
               bool use_caliper_covs,
               NumericVector caliper_covs,
               NumericMatrix caliper_covs_mat,
               bool use_exact,
               IntegerVector exact,
               bool use_antiexact,
               IntegerMatrix antiexact_covs) {

  int k = ii + 1;
  bool found = false, okay;
  int i;

  int n_anti;
  if (use_antiexact) n_anti = antiexact_covs.ncol();

  int n_cal_covs;
  if (use_caliper_covs) n_cal_covs = caliper_covs_mat.ncol();

  while (!found && k <= last_control) {
    if (treat[k] == 1) {
      k++; //if unit is treated, move right
      continue;
    }

    //if unit has already been matched to unit i, skip
    if (r > 0) {
      if (check_in(d_ord[k], row)) {
        k++;
        continue;
      }
    }

    if (!can_be_matched[k]) {
      k++; //if unit is matched, move right
      continue;
    }

    if (use_caliper_dist) {
      if (abs(distance[ii] - distance[k]) > caliper_dist) {
        //if closest is outside caliper, break; none can be found
        break;
      }
    }

    if (use_exact) {
      if (exact[ii] != exact[k]) {
        k++; //if not exact match, move right
        continue;
      }
    }

    if (use_antiexact) {
      i = 0;
      okay = true;
      while (okay && (i < n_anti)) {
        if (antiexact_covs(ii, i) == antiexact_covs(k, i)) {
          okay = false;
        }
        i++;
      }
      if (!okay) {
        k++; //if any antiexact checks failed, move right
        continue;
      }
    }

    if (use_caliper_covs) {
      i = 0;
      okay = true;
      while (okay && (i < n_cal_covs)) {
        if (abs(caliper_covs_mat(ii, i) - caliper_covs_mat(k, i)) > caliper_covs[i]) {
          okay = false;
        }
        i++;
      }
      if (!okay) {
        k++; //if any cov caliper checks failed, move right
        continue;
      }
    }

    found = true;
  }

  if (!found) k = -1;

  return k;
}

int find_left(int ii,
              int first_control,
              IntegerVector treat,
              LogicalVector can_be_matched,
              int r,
              IntegerVector row,
              IntegerVector d_ord,
              NumericVector distance,
              bool use_caliper_dist,
              double caliper_dist,
              bool use_caliper_covs,
              NumericVector caliper_covs,
              NumericMatrix caliper_covs_mat,
              bool use_exact,
              IntegerVector exact,
              bool use_antiexact,
              IntegerMatrix antiexact_covs) {

  int k = ii - 1;
  bool found = false, okay;
  int i;

  int n_anti;
  if (use_antiexact) n_anti = antiexact_covs.ncol();

  int n_cal_covs;
  if (use_caliper_covs) n_cal_covs = caliper_covs_mat.ncol();

  while (!found && k >= first_control) {
    if (treat[k] == 1) {
      k--; //if unit is treated, move left
      continue;
    }

    //if unit has already been matched to unit i, skip
    if (r > 0) {
      if (check_in(d_ord[k], row)) {
        k--;
        continue;
      }
    }

    if (!can_be_matched[k]) {
      k--; //if unit is matched, move left
      continue;
    }

    if (use_caliper_dist) {
      if (abs(distance[ii] - distance[k]) > caliper_dist) {
        //if closest is outside caliper, break
        break;
      }
    }

    if (use_exact) {
      if (exact[ii] != exact[k]) {
        k--; //if not exact match, move left
        continue;
      }
    }

    if (use_antiexact) {
      i = 0;
      okay = true;
      while (okay && (i < n_anti)) {
        if (antiexact_covs(ii, i) == antiexact_covs(k, i)) {
          okay = false;
        }
        i++;
      }
      if (!okay) {
        k--; //if any antiexact checks failed, move left
        continue;
      }
    }

    if (use_caliper_covs) {
      i = 0;
      okay = true;
      while (okay && (i < n_cal_covs)) {
        if (abs(caliper_covs_mat(ii, i) - caliper_covs_mat(k, i)) > caliper_covs[i]) {
          okay = false;
        }
        i++;
      }
      if (!okay) {
        k--; //if any cov caliper checks failed, move left
        continue;
      }
    }

    found = true;
  }

  if (!found) k = -1;

  return k;
}

// [[Rcpp::export]]
IntegerMatrix nn_matchC_vec(const IntegerVector& treat_,
                            const IntegerVector& ord_,
                            const IntegerVector& ratio_,
                            const LogicalVector& discarded_,
                            const int& reuse_max,
                            const NumericVector& distance_,
                            const Nullable<IntegerMatrix>& exact_ = R_NilValue,
                            const Nullable<double>& caliper_dist_ = R_NilValue,
                            const Nullable<NumericVector>& caliper_covs_ = R_NilValue,
                            const Nullable<NumericMatrix>& caliper_covs_mat_ = R_NilValue,
                            const Nullable<IntegerMatrix>& antiexact_covs_ = R_NilValue,
                            const Nullable<IntegerVector>& unit_id_ = R_NilValue,
                            const bool& disl_prog = false) {

  int n = treat_.size();

  CharacterVector lab_ = treat_.names();

  //Use base::sort.list() because faster than Rcpp implementation of order()
  Function o("sort.list");

  IntegerVector d_ord = o(distance_, Named("decreasing") = false);
  d_ord = d_ord - 1;

  IntegerVector treat = treat_[d_ord];
  NumericVector distance = distance_[d_ord];
  CharacterVector lab = lab_[d_ord];
  LogicalVector discarded = discarded_[d_ord];

  IntegerVector ratio_tmp(n);
  ratio_tmp[treat_ == 1] = ratio_;
  IntegerVector ratio = ratio_tmp[d_ord];

  // IntegerVector ord_tmp(n);
  // ord_tmp[treat_ == 1] = ord_;
  // IntegerVector ord = ord_tmp[d_ord];
  // ord = ord - 1;
  IntegerVector ord = ord_ - 1;

  int max_ratio = max(ratio);

  IntegerVector ind = Range(0, n - 1);
  IntegerVector ind0 = ind[treat == 0];
  IntegerVector ind1 = ind[treat == 1];
  int n0 = ind0.size();
  int n1 = ind1.size();
  ind1.names() = lab[ind1];

  IntegerVector t(n);
  IntegerVector t0(n0), t1(n1);
  int i;

  //ind:  1 2 3 4 5 6 7 8
  //ind1:   2 3   5   7

  IntegerMatrix mm(n1, max_ratio);
  mm.fill(NA_INTEGER);
  CharacterVector lab1 = lab[ind1];
  CharacterVector mm_nm = lab_[treat_ == 1];
  rownames(mm) = mm_nm;

  //caliper_dist
  double caliper_dist;
  bool use_caliper_dist = false;
  if (caliper_dist_.isNotNull()) {
    caliper_dist = as<double>(caliper_dist_);
    use_caliper_dist = true;
  }

  //caliper_covs
  NumericVector caliper_covs;
  NumericMatrix caliper_covs_mat;
  bool use_caliper_covs = false;
  if (caliper_covs_.isNotNull()) {
    caliper_covs = as<NumericVector>(caliper_covs_);
    caliper_covs_mat = as<NumericMatrix>(caliper_covs_mat_);
    NumericVector tmp_cc(caliper_covs_mat.nrow());
    for (int i = 0; i < caliper_covs_mat.ncol(); i++) {
      tmp_cc = caliper_covs_mat(_, i);
      tmp_cc = tmp_cc[d_ord];
      caliper_covs_mat(_, i) = tmp_cc;
    }
    use_caliper_covs = true;
  }

  //exact
  bool use_exact = false;
  IntegerVector exact;
  if (exact_.isNotNull()) {
    exact = as<IntegerVector>(exact_)[d_ord];
    use_exact = true;
  }

  //antiexact
  IntegerMatrix antiexact_covs;
  bool use_antiexact = false;
  if (antiexact_covs_.isNotNull()) {
    antiexact_covs = as<IntegerMatrix>(antiexact_covs_);
    NumericVector tmp_ae(antiexact_covs.nrow());
    for (int i = 0; i < antiexact_covs.ncol(); i++) {
      tmp_ae = antiexact_covs(_, i);
      tmp_ae = tmp_ae[d_ord];
      antiexact_covs(_, i) = tmp_ae;
    }
    use_antiexact = true;
  }

  //unit_id
  IntegerVector unit_id;
  bool use_unit_id = false;
  if (unit_id_.isNotNull()) {
    unit_id = as<IntegerVector>(unit_id_)[d_ord];
    use_unit_id = true;
  }

  IntegerVector times_matched = rep(0, n);
  LogicalVector can_be_matched = (!as<LogicalVector>(discarded)) & (treat != 1);

  IntegerVector ind_cbm = ind[can_be_matched];
  int first_control = ind_cbm[0];
  int last_control = ind_cbm[ind_cbm.size() - 1];

  //progress bar
  int prog_length;
  prog_length = max_ratio*n1 + 1;
  Progress p(prog_length, disl_prog);

  int ii, k, j, ck, r, row_to_fill;
  int k_left, k_right;
  double dti;
  String labi;
  IntegerVector mm_row, ck_;

  for (r = 0; r < max_ratio; r++) {
    for (i = 0; i < n1; i++) {
      p.increment();

      row_to_fill = ord[i];

      labi = mm_nm[row_to_fill];

      ii = ind1[labi]; //ii'th unit overall

      if (discarded[ii]) continue;

      if (ratio[ii] < r + 1) continue;

      mm_row = na_omit(mm(row_to_fill, _));

      //find control unit to left and right
      k_left = find_left(ii, first_control,
                         treat,
                         can_be_matched,
                         r, mm_row,
                         d_ord,
                         distance, use_caliper_dist, caliper_dist,
                         use_caliper_covs, caliper_covs, caliper_covs_mat,
                         use_exact, exact,
                         use_antiexact, antiexact_covs);

      k_right = find_right(ii, last_control,
                           treat,
                           can_be_matched,
                           r, mm_row,
                           d_ord,
                           distance, use_caliper_dist, caliper_dist,
                           use_caliper_covs, caliper_covs, caliper_covs_mat,
                           use_exact, exact,
                           use_antiexact,  antiexact_covs);

      if ((k_left >= 0) && (k_right >= 0)) {
        dti = distance[ii];
        if (abs(distance[k_left] - dti) <= abs(distance[k_right] - dti)) {
          k = k_left;
        }
        else {
          k = k_right;
        }
      }
      else if (k_left >= 0) {
        k = k_left;
      }
      else if (k_right >= 0) {
        k = k_right;
      }
      else {
        continue;
      }


      mm( row_to_fill, r ) = d_ord[k] + 1;

      if (use_unit_id) {
        ck_ = ind[unit_id == unit_id[k]];

        for (j = 0; j < ck_.size(); j++) {
          ck = ck_[j];
          times_matched[ck]++;
          if (times_matched[ck] >= reuse_max) {
            can_be_matched[ck] = false;
          }
        }
      }
      else {
        times_matched[k]++;
        if (times_matched[k] >= reuse_max) {
          can_be_matched[k] = false;
        }
      }

      if (any(can_be_matched).is_false()) {
        p.update(prog_length);
        return mm;
      }

      ind_cbm = ind[can_be_matched];
      if (!can_be_matched[first_control]) {
        first_control = ind_cbm[0];
      }
      if (!can_be_matched[last_control]) {
        last_control = ind_cbm[ind_cbm.size() - 1];
      }
    }
  }

  p.update(prog_length);

  return mm;
}
