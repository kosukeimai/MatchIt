// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <Rcpp.h>
#include "internal.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
IntegerMatrix nn_matchC_closest(const NumericMatrix& distance_mat,
                                const IntegerVector& treat,
                                const IntegerVector& ratio,
                                const LogicalVector& discarded,
                                const int& reuse_max,
                                const Nullable<IntegerMatrix>& exact_ = R_NilValue,
                                const Nullable<double>& caliper_dist_ = R_NilValue,
                                const Nullable<NumericVector>& caliper_covs_ = R_NilValue,
                                const Nullable<NumericMatrix>& caliper_covs_mat_ = R_NilValue,
                                const Nullable<IntegerMatrix>& antiexact_covs_ = R_NilValue,
                                const Nullable<IntegerVector>& unit_id_ = R_NilValue,
                                const bool& disl_prog = false)
{

  int r = distance_mat.nrow();
  int c = distance_mat.ncol();

  IntegerMatrix mm(r, max(ratio));
  mm.fill(NA_INTEGER);

  CharacterVector lab = treat.names();

  IntegerVector matched_t = rep(0, r);
  IntegerVector matched_c = rep(0, c);

  // IntegerVector ind = seq(0, treat.size() - 1);
  IntegerVector ind0 = which(treat == 0);
  IntegerVector ind1 = which(treat == 1);

  //caliper_dist
  bool use_caliper_dist = false;
  double caliper_dist;
  if (caliper_dist_.isNotNull()) {
    caliper_dist = as<double>(caliper_dist_);
    use_caliper_dist = true;
  }

  //caliper_covs
  NumericVector caliper_covs;
  NumericMatrix caliper_covs_mat;
  bool use_caliper_covs = false;
  double n_cal_covs;
  if (caliper_covs_.isNotNull()) {
    caliper_covs = as<NumericVector>(caliper_covs_);
    caliper_covs_mat = as<NumericMatrix>(caliper_covs_mat_);
    n_cal_covs = caliper_covs_mat.ncol();
    use_caliper_covs = true;
  }

  //exact
  bool use_exact = false;
  IntegerVector exact;
  if (exact_.isNotNull()) {
    exact = as<IntegerVector>(exact_);
    use_exact = true;
  }

  //antiexact
  IntegerMatrix antiexact_covs;
  bool use_antiexact = false;
  int n_anti;
  if (antiexact_covs_.isNotNull()) {
    antiexact_covs = as<IntegerMatrix>(antiexact_covs_);
    n_anti = antiexact_covs.ncol();
    use_antiexact = true;
  }

  //unit_id
  IntegerVector unit_id, ck_;
  bool use_unit_id = false;
  if (unit_id_.isNotNull()) {
    unit_id = as<IntegerVector>(unit_id_);
    use_unit_id = true;
  }

  //progress bar
  int prog_length;
  prog_length = sum(ratio) + 1;
  Progress p(prog_length, disl_prog);
  p.increment();

  Function o("order");

  IntegerVector d_ord = o(distance_mat);
  d_ord = d_ord - 1; //Because R uses 1-indexing

  int rj, cj, dj, i, ind0i, ind1i;
  bool okay;

  for (int j = 0; j < d_ord.size(); j++) {

    dj = d_ord[j];

    // If distance is greater tha distance caliper, stop the whole thing because
    // no remaining distance will be smaller
    if (use_caliper_dist) {
      if (distance_mat[dj] > caliper_dist) break;
    }

    // Get row and column index of potential pair
    rj = dj % r;
    cj = dj / r;

    // Get sample indices of members of potential pair
    ind1i = ind1[rj];
    ind0i = ind0[cj];

    // If either member is discarded, move on
    if (discarded[ind1i]) continue;
    if (discarded[ind0i]) continue;

    // If either member has been matched enough times, move on
    if (matched_t[rj] >= ratio[rj]) continue;
    if (matched_c[cj] >= reuse_max) continue;

    // Exact matching criterion
    if (use_exact) {
      if (exact[ind1i] != exact[ind0i]) {
        continue;
      }
    }

    // Covariate caliper criterion
    if (use_caliper_covs) {
      i = 0;
      okay = true;
      while (okay && (i < n_cal_covs)) {
        if (std::abs(caliper_covs_mat(ind1i, i) - caliper_covs_mat(ind0i, i)) > caliper_covs[i]) {
          okay = false;
        }
        i++;
      }
      if (!okay) continue;
    }

    // Antiexact criterion
    if (use_antiexact) {
      i = 0;
      okay = true;
      while (okay && (i < n_anti)) {
        if (antiexact_covs(ind1i, i) == antiexact_covs(ind0i, i)) {
          okay = false;
        }
        i++;
      }
      if (!okay) continue;
    }

    // If all criteria above are satisfied, potential pair becomes a pair!

    // If unit_id used, increase match count of all units with that ID
    if (use_unit_id) {
      ck_ = which(as<IntegerVector>(unit_id[ind1]) == unit_id[ind1i]);

      for (i = 0; i < ck_.size(); i++) {
        matched_t[ck_[i]]++;
      }

      ck_ = which(as<IntegerVector>(unit_id[ind0]) == unit_id[ind0i]);

      for (i = 0; i < ck_.size(); i++) {
        matched_c[ck_[i]]++;
      }
    }
    else {
      matched_t[rj]++;
      matched_c[cj]++;
    }

    mm(rj, matched_t[rj] - 1) = ind0i;

    p.increment();

    if (matched_t[rj] >= ratio[rj]) {
      if (all(matched_t >= ratio).is_true()) break;
    }
    if (matched_c[cj] >= reuse_max) {
      if (all(matched_c >= reuse_max).is_true()) break;
    }
  }

  p.update(prog_length);

  mm = mm + 1;
  rownames(mm) = lab[treat == 1];

  return mm;
}
