// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include "eta_progress_bar.h"
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
                                const bool& disl_prog = false) {

  IntegerVector unique_treat = {0, 1};
  int g = unique_treat.size();
  int focal = 1;

  R_xlen_t n = treat.size();
  IntegerVector ind = Range(0, n - 1);

  R_xlen_t i;
  int gi;
  IntegerVector indt(n);
  IntegerVector indt_sep(g + 1);
  IntegerVector indt_tmp;
  IntegerVector nt(g);
  IntegerVector ind_match(n);
  ind_match.fill(NA_INTEGER);

  IntegerVector times_matched(n);
  times_matched.fill(0);
  LogicalVector eligible = !discarded;

  for (gi = 0; gi < g; gi++) {
    nt[gi] = sum(treat == gi);
  }

  int nf = nt[focal];

  indt_sep[0] = 0;

  for (gi = 0; gi < g; gi++) {
    indt_sep[gi + 1] = indt_sep[gi] + nt[gi];

    indt_tmp = ind[treat == gi];

    for (i = 0; i < nt[gi]; i++) {
      indt[indt_sep[gi] + i] = indt_tmp[i];
      ind_match[indt_tmp[i]] = i;
    }
  }

  IntegerVector ind_focal = indt[Range(indt_sep[focal], indt_sep[focal + 1] - 1)];

  IntegerVector times_matched_allowed(n);
  times_matched_allowed.fill(reuse_max);
  times_matched_allowed[ind_focal] = ratio;

  IntegerVector n_eligible(g);
  for (i = 0; i < n; i++) {
    if (eligible[i]) {
      n_eligible[treat[i]]++;
    }
  }

  int max_ratio = max(ratio);

  // Output matrix with sample indices of control units
  IntegerMatrix mm(nf, max_ratio);
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
    caliper_dist = max_finite(distance_mat) + .1;
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

  //progress bar
  int prog_length;
  prog_length = sum(ratio) + 1;
  ETAProgressBar pb;
  Progress p(prog_length, disl_prog, pb);

  Function o("order");

  IntegerVector d_ord = o(distance_mat);
  d_ord = d_ord - 1; //Because R uses 1-indexing

  gi = 0;

  R_xlen_t r = distance_mat.nrow();

  int rj, cj, c_id_i, t_id_i;
  int counter = -1;

  for (R_xlen_t dj : d_ord) {
    counter++;
    if (counter % 200 == 0) Rcpp::checkUserInterrupt();

    if (min(n_eligible) <= 0) {
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
    t_id_i = ind_focal[rj];

    // If either member is discarded, move on
    if (!eligible[t_id_i]) {
      continue;
    }

    c_id_i = indt[indt_sep[gi] + cj];

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
    ck_ = {t_id_i, c_id_i};

    if (use_unit_id) {
      ck_ = which(!is_na(match(unit_id, as<IntegerVector>(unit_id[ck_]))));
    }

    for (int ck : ck_) {

      if (!eligible[ck]) {
        continue;
      }

      times_matched[ck]++;
      if (times_matched[ck] >= times_matched_allowed[ck]) {
        eligible[ck] = false;
        n_eligible[treat[ck]]--;
      }
    }

    p.increment();
  }

  p.update(prog_length);

  mm = mm + 1;
  rownames(mm) = lab[ind_focal];

  return mm;
}
