// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include "eta_progress_bar.h"
#include <Rcpp.h>
#include "internal.h"
#include <cmath>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
IntegerMatrix nn_matchC(const IntegerVector& treat_,
                        const IntegerVector& ord,
                        const IntegerVector& ratio,
                        const LogicalVector& discarded,
                        const int& reuse_max,
                        const int& focal_,
                        const Nullable<NumericVector>& distance_ = R_NilValue,
                        const Nullable<NumericMatrix>& distance_mat_ = R_NilValue,
                        const Nullable<IntegerVector>& exact_ = R_NilValue,
                        const Nullable<double>& caliper_dist_ = R_NilValue,
                        const Nullable<NumericVector>& caliper_covs_ = R_NilValue,
                        const Nullable<NumericMatrix>& caliper_covs_mat_ = R_NilValue,
                        const Nullable<NumericMatrix>& mah_covs_ = R_NilValue,
                        const Nullable<IntegerMatrix>& antiexact_covs_ = R_NilValue,
                        const Nullable<IntegerVector>& unit_id_ = R_NilValue,
                        const bool& disl_prog = false) {

  IntegerVector unique_treat = unique(treat_);
  std::sort(unique_treat.begin(), unique_treat.end());
  int g = unique_treat.size();
  IntegerVector treat = match(treat_, unique_treat) - 1;
  int focal;
  for (focal = 0; focal < g; focal++) {
    if (unique_treat[focal] == focal_) {
      break;
    }
  }

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

  IntegerVector g_c = Range(0, g - 1);
  g_c = g_c[g_c != focal];

  for (gi = 0; gi < g; gi++) {
    nt[gi] = sum(treat == gi);
  }

  int nf = nt[focal];

  int max_nc = max(as<IntegerVector>(nt[g_c]));

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

  IntegerVector n_eligible(unique_treat.size());
  for (i = 0; i < n; i++) {
    if (eligible[i]) {
      n_eligible[treat[i]]++;
    }
  }

  int max_ratio = max(ratio);

  // Output matrix with sample indices of control units
  IntegerMatrix mm(nf, max_ratio);
  mm.fill(NA_INTEGER);
  CharacterVector lab = treat_.names();

  int min_ind, t_rat;

  //exact
  bool use_exact = false;
  IntegerVector exact;
  if (exact_.isNotNull()) {
    exact = as<IntegerVector>(exact_);
    use_exact = true;
  }

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
  int ncc = 0;
  if (caliper_covs_.isNotNull()) {
    caliper_covs = as<NumericVector>(caliper_covs_);
    caliper_covs_mat = as<NumericMatrix>(caliper_covs_mat_);
    ncc = caliper_covs_mat.ncol();
  }

  //dsit_mat and mah_covs
  bool use_dist_mat = false;
  bool use_mah_covs = false;
  NumericMatrix distance_mat, mah_covs;
  if (mah_covs_.isNotNull()) {
    mah_covs = as<NumericMatrix>(mah_covs_);
    use_mah_covs = true;
  }
  else if (distance_mat_.isNotNull()) {
    distance_mat = as<NumericMatrix>(distance_mat_);
    use_dist_mat = true;
  }

  //distance
  NumericVector distance;
  if (distance_.isNotNull()) {
    distance = distance_;
  }

  //anitexact
  IntegerMatrix antiexact_covs;
  int aenc = 0;
  if (antiexact_covs_.isNotNull()) {
    antiexact_covs = as<IntegerMatrix>(antiexact_covs_);
    aenc = antiexact_covs.ncol();
  }

  //reuse_max
  bool use_reuse_max = (reuse_max < nf);

  //unit_id
  IntegerVector unit_id;
  IntegerVector matched_unit_ids;
  bool use_unit_id = false;
  if (unit_id_.isNotNull()) {
    unit_id = as<IntegerVector>(unit_id_);
    use_unit_id = true;
    use_reuse_max = true;
    matched_unit_ids = rep(NA_INTEGER, max_nc);
  }

  IntegerVector c_eligible(max_nc);
  NumericVector match_distance(max_nc);
  IntegerVector matches_i(1 + max_ratio * (g - 1));
  int k_total;

  //progress bar
  int prog_length;
  if (use_reuse_max) prog_length = sum(ratio) + 1;
  else prog_length = nf + 1;
  ETAProgressBar pb;
  Progress p(prog_length, disl_prog, pb);

  //Counters
  int r, t_id_t_i, t_id_i, c_id_i, c, k;
  double ps_diff, dist;
  IntegerVector ck_, top_r_matches;
  bool ps_diff_calculated;

  int counter = -1;

  //Matching
  if (use_reuse_max) {
    for (r = 1; r <= max_ratio; r++) {
      for (i = 0; i < nf && max(as<IntegerVector>(n_eligible[g_c])) > 0; i++) {

        counter++;
        if (counter % 200 == 0) Rcpp::checkUserInterrupt();

        t_id_t_i = ord[i] - 1; // index among treated
        t_id_i = ind_focal[t_id_t_i]; // index among sample

        if (r > times_matched_allowed[t_id_i]) {
          continue;
        }

        p.increment();

        if (!eligible[t_id_i]) {
          continue;
        }

        k_total = 0;

        for (int gi : g_c) {
          k = 0;

          if (n_eligible[gi] > 0) {
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

              ps_diff_calculated = false;

              if (use_caliper_dist) {
                if (use_dist_mat) {
                  ps_diff = distance_mat(t_id_t_i, ind_match[c_id_i]);
                }
                else {
                  ps_diff = std::abs(distance[c_id_i] - distance[t_id_i]);
                }

                if (ps_diff > caliper_dist) {
                  continue;
                }

                ps_diff_calculated = true;
              }

              if (!caliper_covs_okay(ncc, t_id_i, c_id_i, caliper_covs_mat, caliper_covs)) {
                continue;
              }

              //Compute distances among eligible
              if (use_mah_covs) {
                dist = sum(pow(mah_covs.row(t_id_i) - mah_covs.row(c_id_i), 2.0));
              }
              else if (ps_diff_calculated) {
                dist = ps_diff;
              }
              else if (use_dist_mat) {
                dist = distance_mat(t_id_t_i, ind_match[c_id_i]);
              }
              else {
                dist = std::abs(distance[c_id_i] - distance[t_id_i]);
              }

              if (!std::isfinite(dist)) {
                continue;
              }

              c_eligible[k] = c_id_i;
              match_distance[k] = dist;
              k++;
            }
          }

          //If no matches...
          if (k == 0) {
            //If round 1, focal has no possible matches
            if (r == 1) {
              k_total = 0;
              break;
            }
            continue;
          }

          //Find minimum distance and assign
          min_ind = 0;
          for (c = 1; c < k; c++) {
            if (match_distance[c] < match_distance[min_ind]) {
              min_ind = c;
            }
          }

          matches_i[k_total] = c_eligible[min_ind];
          k_total++;
        }

        if (k_total == 0) {
          eligible[t_id_i] = false;
          n_eligible[focal]--;
          continue;
        }

        //Assign to match matrix
        for (c = 0; c < k_total; c++) {
          mm(t_id_t_i, sum(!is_na(mm(t_id_t_i, _)))) = matches_i[c];
        }

        matches_i[k_total] = t_id_i;

        ck_ = matches_i[Range(0, k_total)];

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
      }
    }
  }
  else {
    for (i = 0; i < nf; i++) {

      counter++;
      if (counter % 500 == 0) Rcpp::checkUserInterrupt();

      t_id_t_i = ord[i] - 1; // index among treated
      t_id_i = ind_focal[t_id_t_i]; // index among sample

      p.increment();

      if (!eligible[t_id_i]) {
        continue;
      }

      t_rat = ratio[t_id_t_i];

      k_total = 0;

      for (int gi : g_c) {
        k = 0;

        if (n_eligible[gi] > 0) {
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

            ps_diff_calculated = false;

            if (use_caliper_dist) {
              if (use_dist_mat) {
                ps_diff = distance_mat(t_id_t_i, ind_match[c_id_i]);
              }
              else {
                ps_diff = std::abs(distance[c_id_i] - distance[t_id_i]);
              }

              if (ps_diff > caliper_dist) {
                continue;
              }

              ps_diff_calculated = true;
            }

            if (!caliper_covs_okay(ncc, t_id_i, c_id_i, caliper_covs_mat, caliper_covs)) {
              continue;
            }

            //Compute distances among eligible
            if (use_mah_covs) {
              dist = sqrt(sum(pow(mah_covs.row(t_id_i) - mah_covs.row(c_id_i), 2.0)));
            }
            else if (ps_diff_calculated) {
              dist = ps_diff;
            }
            else if (use_dist_mat) {
              dist = distance_mat(t_id_t_i, ind_match[c_id_i]);
            }
            else {
              dist = std::abs(distance[c_id_i] - distance[t_id_i]);
            }

            if (!std::isfinite(dist)) {
              continue;
            }

            c_eligible[k] = c_id_i;
            match_distance[k] = dist;
            k++;
          }
        }

        //If no matches...
        if (k == 0) {
          k_total = 0;
          break;
        }

        //If replace and few eligible controls, assign all and move on

        if (k < t_rat) {
          t_rat = k;
        }

        //Sort distances and assign
        top_r_matches = Range(0, k - 1);

        std::partial_sort(top_r_matches.begin(), top_r_matches.begin() + t_rat, top_r_matches.end(),
                          [&match_distance](int a, int b) {return match_distance[a] < match_distance[b];});

        for (c = 0; c < t_rat; c++) {
          matches_i[k_total] = c_eligible[top_r_matches[c]];
          k_total++;
        }
      }

      if (k_total == 0) {
        continue;
      }

      //Assign to match matrix
      for (c = 0; c < k_total; c++) {
        mm(t_id_t_i, sum(!is_na(mm(t_id_t_i, _)))) = matches_i[c];
      }
    }
  }

  p.update(prog_length);

  mm = mm + 1; // + 1 because C indexing starts at 0 but mm is sent to R
  rownames(mm) = lab[ind_focal];

  return mm;
}
