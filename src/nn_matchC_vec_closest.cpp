// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <Rcpp.h>
#include "internal.h"
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
IntegerMatrix nn_matchC_vec_closest(const IntegerVector& treat,
                                    const IntegerVector& ratio,
                                    const LogicalVector& discarded,
                                    const int& reuse_max,
                                    const NumericVector& distance,
                                    const Nullable<IntegerMatrix>& exact_ = R_NilValue,
                                    const Nullable<double>& caliper_dist_ = R_NilValue,
                                    const Nullable<NumericVector>& caliper_covs_ = R_NilValue,
                                    const Nullable<NumericMatrix>& caliper_covs_mat_ = R_NilValue,
                                    const Nullable<IntegerMatrix>& antiexact_covs_ = R_NilValue,
                                    const Nullable<IntegerVector>& unit_id_ = R_NilValue,
                                    const bool& disl_prog = false) {

  int n = treat.size();

  //Use base::order() because faster than Rcpp implementation of order()
  Function o("order");

  IntegerVector ind_d_ord = o(distance, Named("decreasing") = false);
  ind_d_ord = ind_d_ord - 1; //location of each unit after sorting

  IntegerVector ind = Range(0, n - 1);
  IntegerVector ind1 = ind[treat == 1];

  int n1 = ind1.size();

  int i, j;

  IntegerVector match_d_ord = match(ind, ind_d_ord) - 1;

  int max_ratio = max(ratio);

  //ind:  1 2 3 4 5 6 7 8
  //ind1:   2 3   5   7

  IntegerMatrix mm(n1, max_ratio);
  mm.fill(NA_INTEGER);
  CharacterVector lab = treat.names();

  //exact
  bool use_exact = false;
  IntegerVector exact;
  if (exact_.isNotNull()) {
    exact = as<IntegerVector>(exact_);
    use_exact = true;
  }

  //caliper_dist
  double caliper_dist;
  if (caliper_dist_.isNotNull()) {
    caliper_dist = as<double>(caliper_dist_);
  }
  else {
    caliper_dist = max_finite(distance) - min_finite(distance) + 1;
  }

  //caliper_covs
  NumericVector caliper_covs;
  NumericMatrix caliper_covs_mat;
  int ncc = 0;
  if (caliper_covs_.isNotNull()) {
    caliper_covs = as<NumericVector>(caliper_covs_);
    caliper_covs_mat = as<NumericMatrix>(caliper_covs_mat_);
    ncc = caliper_covs_mat.ncol();
  }

  //antiexact
  IntegerMatrix antiexact_covs;
  int aenc = 0;
  if (antiexact_covs_.isNotNull()) {
    antiexact_covs = as<IntegerMatrix>(antiexact_covs_);
    aenc = antiexact_covs.ncol();
  }

  //unit_id
  IntegerVector unit_id;
  bool use_unit_id = false;
  if (unit_id_.isNotNull()) {
    unit_id = as<IntegerVector>(unit_id_);
    use_unit_id = true;
  }

  IntegerVector ind1_match = rep(NA_INTEGER, n);
  for (i = 0; i < n1; i++) {
    ind1_match[ind1[i]] = i;
  }

  //storing closeness
  IntegerVector t_id(2 * n1);
  IntegerVector c_id(2 * n1);
  NumericVector dist = rep(R_PosInf, 2 * n1);
  for (i = 0; i < n1; i++) {
    t_id[i] = ind1[i];
    t_id[i + n1] = ind1[i];
    c_id[i] = -1;
    c_id[i + n1] = -2;
  }

  IntegerVector times_matched = rep(0, n);
  LogicalVector eligible = rep(true, n);
  eligible[discarded] = false;
  IntegerVector times_matched_allowed = rep(reuse_max, n);
  times_matched_allowed[ind1] = ratio;

  IntegerVector times_skipped = rep(0, n);

  //progress bar
  int prog_length = sum(ratio) + 1;
  Progress p(prog_length, disl_prog);
  p.increment();

  IntegerVector ck_;

  int t_id_i, c_id_i, t_id_t_i, k;

  int first_control = 0;
  int last_control = n - 1;

  while (!eligible[ind_d_ord[first_control]] || treat[ind_d_ord[first_control]] == 1) {
    first_control++;
  }
  while (!eligible[ind_d_ord[last_control]] || treat[ind_d_ord[last_control]] == 1) {
    last_control--;
  }

  //Find left and right matches for each treated unit

  for (i = 0; i < (2 * n1); i++) {
    t_id_i = t_id[i];

    if (!eligible[t_id_i]) {
      continue;
    }

    c_id_i = c_id[i];

    k = find_lr(c_id_i,
                t_id_i,
                ind_d_ord,
                match_d_ord,
                treat,
                distance,
                eligible,
                0,
                ncc,
                caliper_covs_mat,
                caliper_covs,
                caliper_dist,
                use_exact,
                exact,
                aenc,
                antiexact_covs,
                first_control,
                last_control);

    if (k < 0) {
      times_skipped[t_id_i]++;
      if (times_skipped[t_id_i] == 2) {
        eligible[t_id_i] = false;
      }
      continue;
    }

    c_id[i] = k;
    dist[i] = std::abs(distance[t_id_i] - distance[k]);
  }

  int n_eligible0 = 0;
  int n_eligible1 = 0;
  for (i = 0; i < n; i++) {
    if (eligible[i]) {
      if (treat[i] == 0) {
        n_eligible0++;
      }
      else {
        n_eligible1++;
      }
    }
  }

  //Order the list
  IntegerVector heap_ord = o(dist);
  heap_ord = heap_ord - 1;

  //Go down the list; update as needed
  int hi;
  int counter = -1;

  for (i = 0; (i < 2 * n1) && n_eligible1 > 0 && n_eligible0 > 0; i++) {
    counter++;
    if (counter % 500 == 0) Rcpp::checkUserInterrupt();

    hi = heap_ord[i];

    if (dist[hi] >= caliper_dist) {
      break;
    }

    t_id_i = t_id[hi];

    if (!eligible[t_id_i]) {
      continue;
    }

    c_id_i = c_id[hi];

    if (!eligible[c_id_i]) {
      // If control isn't eligible, find new control and try again

      while (!eligible[ind_d_ord[first_control]] || treat[ind_d_ord[first_control]] == 1) {
        first_control++;
      }
      while (!eligible[ind_d_ord[last_control]] || treat[ind_d_ord[last_control]] == 1) {
        last_control--;
      }

      k = find_lr(c_id_i,
                  t_id_i,
                  ind_d_ord,
                  match_d_ord,
                  treat,
                  distance,
                  eligible,
                  0,
                  ncc,
                  caliper_covs_mat,
                  caliper_covs,
                  caliper_dist,
                  use_exact,
                  exact,
                  aenc,
                  antiexact_covs,
                  first_control,
                  last_control);

      //If no new control found, mark treated unit as ineligible and continue
      if (k < 0) {
        times_skipped[t_id_i]++;
        if (times_skipped[t_id_i] == 2) {
          eligible[t_id_i] = false;
          n_eligible1--;
        }
        continue;
      }

      c_id[hi] = k;
      dist[hi] = std::abs(distance[t_id_i] - distance[k]);

      //Find new position of pair in heap
      for (j = i; j < (2 * n1) - 1; j++) {
        if (dist[heap_ord[j]] < dist[heap_ord[j + 1]]) {
          break;
        }

        swap_pos(heap_ord, j, j + 1);
      }

      i--;
    }
    else {
      t_id_t_i = ind1_match[t_id_i];

      mm(t_id_t_i, sum(!is_na(mm(t_id_t_i, _)))) = c_id_i;

      if (use_unit_id) {
        ck_ = ind[unit_id == unit_id[t_id_i] | unit_id == unit_id[c_id_i]];
      }
      else {
        ck_ = {c_id_i, t_id_i};
      }

      for (int ck : ck_) {

        if (!eligible[ck]) {
          continue;
        }

        times_matched[ck]++;
        if (times_matched[ck] >= times_matched_allowed[ck]) {
          eligible[ck] = false;
          if (treat[ck] == 1) {
            n_eligible1--;
          }
          else {
            n_eligible0--;
          }
        }
      }

      p.increment();
    }
  }

  p.update(prog_length);

  mm = mm + 1;
  rownames(mm) = lab[ind1];

  return mm;
}