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

  int n = treat.size();

  IntegerMatrix mm(r, max(ratio));
  mm.fill(NA_INTEGER);

  CharacterVector lab = treat.names();


  IntegerVector ind = Range(0, n - 1);
  IntegerVector ind0 = ind[treat == 0];
  IntegerVector ind1 = ind[treat == 1];

  //caliper_dist
  double caliper_dist;
  if (caliper_dist_.isNotNull()) {
    caliper_dist = as<double>(caliper_dist_);
  }
  else {
    caliper_dist = max_finite(distance_mat) + 1;
  }

  //caliper_covs
  NumericVector caliper_covs;
  NumericMatrix caliper_covs_mat;
  int ncc;
  if (caliper_covs_.isNotNull()) {
    caliper_covs = as<NumericVector>(caliper_covs_);
    caliper_covs_mat = as<NumericMatrix>(caliper_covs_mat_);
    ncc = caliper_covs_mat.ncol();
  }
  else {
    ncc = 0;
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
  int aenc;
  if (antiexact_covs_.isNotNull()) {
    antiexact_covs = as<IntegerMatrix>(antiexact_covs_);
    aenc = antiexact_covs.ncol();
  }
  else {
    aenc = 0;
  }

  //unit_id
  IntegerVector unit_id, ck_;
  bool use_unit_id = false;
  if (unit_id_.isNotNull()) {
    unit_id = as<IntegerVector>(unit_id_);
    use_unit_id = true;
  }

  IntegerVector times_matched = rep(0, n);
  LogicalVector eligible = rep(true, n);
  eligible[discarded] = false;
  IntegerVector times_matched_allowed = rep(reuse_max, n);
  times_matched_allowed[ind1] = ratio;

  int n_eligible0 = sum(as<LogicalVector>(eligible[treat == 0]));
  int n_eligible1 = sum(as<LogicalVector>(eligible[treat == 1]));

  //progress bar
  int prog_length;
  prog_length = sum(ratio) + 1;
  Progress p(prog_length, disl_prog);
  p.increment();

  Function o("order");

  IntegerVector d_ord = o(distance_mat);
  d_ord = d_ord - 1; //Because R uses 1-indexing

  int rj, cj, c_id_i, t_id_i;

  for (int dj : d_ord) {

    if (n_eligible1 <= 0) {
      break;
    }

    if (n_eligible0 <= 0) {
      break;
    }

    // If distance is greater than distance caliper, stop the whole thing because
    // no remaining distance will be smaller
    if (distance_mat[dj] > caliper_dist) {
      break;
    }

    // Get row and column index of potential pair
    rj = dj % r;
    cj = dj / r;

    // Get sample indices of members of potential pair
    t_id_i = ind1[rj];
    c_id_i = ind0[cj];

    // If either member is discarded, move on
    if (!eligible[t_id_i]) {
      continue;
    }

    if (!eligible[c_id_i]) {
      continue;
    }

    // Exact matching criterion
    if (!exact_okay(use_exact, t_id_i, c_id_i, exact)) {
      continue;
    }

    // Antiexact criterion
    if (!antiexact_okay(aenc, t_id_i, c_id_i, antiexact_covs)) {
      continue;
    }

    // Covariate caliper criterion
    if (!caliper_covs_okay(ncc, t_id_i, c_id_i, caliper_covs_mat, caliper_covs)) {
      continue;
    }

    // If all criteria above are satisfied, potential pair becomes a pair!
    mm(rj, sum(!is_na(mm(rj, _)))) = c_id_i;

    // If unit_id used, increase match count of all units with that ID
    if (use_unit_id) {
      ck_ = ind[unit_id == unit_id[t_id_i] | unit_id == unit_id[c_id_i]];
    }
    else {
      ck_ = {t_id_i, c_id_i};
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

  p.update(prog_length);

  mm = mm + 1;
  rownames(mm) = lab[treat == 1];

  return mm;
}
