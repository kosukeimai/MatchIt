// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include "eta_progress_bar.h"
#include <Rcpp.h>
#include "internal.h"
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
IntegerMatrix nn_matchC_mahcovs_closest(const IntegerVector& treat,
                                        const IntegerVector& ratio,
                                        const LogicalVector& discarded,
                                        const int& reuse_max,
                                        const NumericMatrix& mah_covs,
                                        const Nullable<NumericVector>& distance_ = R_NilValue,
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

  //distance
  bool use_caliper_dist = false;
  double caliper_dist, ps_diff;
  NumericVector distance;
  if (distance_.isNotNull()) {
    distance = distance_;

    //caliper_dist
    if (caliper_dist_.isNotNull()) {
      caliper_dist = as<double>(caliper_dist_);
      use_caliper_dist = true;
    }
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

  //storing closeness
  IntegerVector t_id = ind_focal;
  IntegerVector c_id = rep(-1, nf);
  NumericVector dist = rep(R_PosInf, nf);

  //progress bar
  int prog_length = nf + sum(ratio) + 1;
  ETAProgressBar pb;
  Progress p(prog_length, disl_prog, pb);

  gi = 0;

  IntegerVector ck_, t_inds;
  IntegerVector c_eligible(nt[gi]);
  NumericVector match_distance(nt[gi]);

  R_xlen_t c;
  int c_id_i, t_id_t_i;
  int t_id_i = -1;
  double dist_c;
  bool any_match_found;

  int counter = -1;
  int r = 1;

  //Find closest control unit to each treated unit
  for (i = 0; i < nf; i++) {
    counter++;
    if (counter % 200 == 0) Rcpp::checkUserInterrupt();

    p.increment();

    t_id_i = ind_focal[i];

    if (!eligible[t_id_i]) {
      continue;
    }

    t_inds = which(t_id == t_id_i);

    any_match_found = false;

    for (c = indt_sep[gi]; c < indt_sep[gi + 1]; c++) {
      c_id_i = indt[c];

      if (!eligible[c_id_i]) {
        continue;
      }

      if (!exact_okay(use_exact, t_id_i, c_id_i, exact)) {
        continue;
      }

      if (!antiexact_okay(aenc, t_id_i, c_id_i, antiexact_covs)) {
        continue;
      }

      if (use_caliper_dist) {
        ps_diff = std::abs(distance[c_id_i] - distance[t_id_i]);

        if (ps_diff > caliper_dist) {
          continue;
        }
      }

      if (!caliper_covs_okay(ncc, t_id_i, c_id_i, caliper_covs_mat, caliper_covs)) {
        continue;
      }

      //Compute distances among eligible
      dist_c = sqrt(sum(pow(mah_covs.row(t_id_i) - mah_covs.row(c_id_i), 2.0)));

      if (!std::isfinite(dist_c)) {
        continue;
      }

      if (any_match_found) {
        if (dist_c < dist[i]) {
          c_id[i] = c_id_i;
          dist[i] = dist_c;
        }
      }
      else {
        c_id[i] = c_id_i;
        dist[i] = dist_c;
        any_match_found = true;
      }
    }

    if (!any_match_found) {
      eligible[t_id_i] = false;
      n_eligible[focal]--;
    }
  }

  //Order the list
  // Use base::order() because faster than Rcpp implementation of order()
  Function o("order");
  IntegerVector heap_ord = o(dist);
  heap_ord = heap_ord - 1;

  //Go down the list; update as needed
  R_xlen_t hi;
  bool find_new;

  i = 0;
  while (min(n_eligible) > 0 && i < nf) {
    counter++;
    if (counter % 200 == 0) Rcpp::checkUserInterrupt();

    hi = heap_ord[i];

    t_id_i = t_id[hi];

    if (!eligible[t_id_i]) {
      i++;
      continue;
    }

    r = times_matched[t_id_i] + 1;

    t_id_t_i = ind_match[t_id_i];

    c_id_i = c_id[hi];

    find_new = false;
    if (!eligible[c_id_i]) {
      find_new = true;
    }
    else if (!mm_okay(r, c_id_i, mm.row(t_id_t_i))) {
      find_new = true;
    }

    if (find_new) {
      // If control isn't eligible, find new control and try again
      any_match_found = false;

      for (c = indt_sep[gi]; c < indt_sep[gi + 1]; c++) {
        c_id_i = indt[c];

        if (!eligible[c_id_i]) {
          continue;
        }

        //Prevent control units being matched to same treated unit again
        if (!mm_okay(r, c_id_i, mm.row(t_id_t_i))) {
          continue;
        }

        if (!exact_okay(use_exact, t_id_i, c_id_i, exact)) {
          continue;
        }

        if (!antiexact_okay(aenc, t_id_i, c_id_i, antiexact_covs)) {
          continue;
        }

        if (use_caliper_dist) {
          ps_diff = std::abs(distance[c_id_i] - distance[t_id_i]);

          if (ps_diff > caliper_dist) {
            continue;
          }
        }

        if (!caliper_covs_okay(ncc, t_id_i, c_id_i, caliper_covs_mat, caliper_covs)) {
          continue;
        }

        //Compute distances among eligible
        dist_c = sum(pow(mah_covs.row(t_id_i) - mah_covs.row(c_id_i), 2.0));

        if (!std::isfinite(dist_c)) {
          continue;
        }

        if (any_match_found) {
          if (dist_c < dist[hi]) {
            c_id[hi] = c_id_i;
            dist[hi] = dist_c;
          }
        }
        else {
          c_id[hi] = c_id_i;
          dist[hi] = dist_c;
          any_match_found = true;
        }
      }

      //If no matches...
      if (!any_match_found) {
        eligible[t_id_i] = false;
        n_eligible[focal]--;
        continue;
      }

      //Find new position of pair in heap
      for (c = i; c < nf - 1; c++) {
        if (dist[heap_ord[c]] < dist[heap_ord[c + 1]]) {
          break;
        }

        swap_pos(heap_ord, c, c + 1);
      }

      continue;
    }

    mm(t_id_t_i, sum(!is_na(mm(t_id_t_i, _)))) = c_id_i;

    ck_ = {c_id_i, t_id_i};

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
