// [[Rcpp::depends(RcppProgress)]]
#include "eta_progress_bar.h"
#include "internal.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
IntegerMatrix nn_matchC_distmat(const IntegerVector& treat_,
                                const IntegerVector& ord,
                                const IntegerVector& ratio,
                                const LogicalVector& discarded,
                                const int& reuse_max,
                                const int& focal_,
                                const NumericMatrix& distance_mat,
                                const Nullable<IntegerMatrix>& exact_ = R_NilValue,
                                const Nullable<double>& caliper_dist_ = R_NilValue,
                                const Nullable<NumericVector>& caliper_covs_ = R_NilValue,
                                const Nullable<NumericMatrix>& caliper_covs_mat_ = R_NilValue,
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

  LogicalVector eligible = !discarded;

  IntegerVector g_c = Range(0, g - 1);
  g_c = g_c[g_c != focal];

  IntegerVector n_eligible(g);
  for (i = 0; i < n; i++) {
    nt[treat[i]]++;

    if (eligible[i]) {
      n_eligible[treat[i]]++;
    }
  }

  int nf = nt[focal];

  indt_sep[0] = 0;

  for (gi = 0; gi < g; gi++) {
    indt_sep[gi + 1] = indt_sep[gi] + nt[gi];

    indt_tmp = ind[treat == gi];

    for (i = 0; i < nt[gi]; i++) {
      indt[indt_sep[gi] + i] = indt_tmp[i];
    }
  }

  IntegerVector ind_focal = indt[Range(indt_sep[focal], indt_sep[focal + 1] - 1)];

  std::vector<int> times_matched(n, 0);

  std::vector<int> times_matched_allowed(n, reuse_max);
  for (i = 0; i < nf; i++) {
    times_matched_allowed[ind_focal[i]] = ratio[i];
  }

  int max_ratio = max(ratio);

  IntegerVector ind_non_focal = which(treat != focal);

  for (i = 0; i < n - nf; i++) {
    ind_match[ind_non_focal[i]] = i;
  }

  for (i = 0; i < nf; i++) {
    ind_match[ind_focal[i]] = i;
  }

  // Output matrix with sample indices of control units
  IntegerMatrix mm(nf, max_ratio);
  mm.fill(NA_INTEGER);
  CharacterVector lab = treat_.names();

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

  //reuse_max
  bool use_reuse_max = (reuse_max < nf);

  //unit_id
  IntegerVector unit_id;
  bool use_unit_id = false;
  if (unit_id_.isNotNull()) {
    unit_id = as<IntegerVector>(unit_id_);
    use_unit_id = true;
    use_reuse_max = true;
  }

  IntegerVector matches_i(1 + max_ratio * (g - 1));
  int k_total;

  //progress bar
  int prog_length;
  if (use_reuse_max) prog_length = sum(ratio) + 1;
  else prog_length = nf + 1;
  ETAProgressBar pb;
  Progress p(prog_length, disl_prog, pb);

  R_xlen_t c;
  int r, t_id_i;
  IntegerVector ck_;
  std::vector<int> k(max_ratio);

  int counter = 0;

  if (use_reuse_max) {
    IntegerVector ord_r(nf);

    for (r = 1; r <= max_ratio; r++) {
      ord_r = ord[as<IntegerVector>(ratio[ord - 1]) >= r];
      ord_r = ord_r - 1;

      for (int t_id_t_i : ord_r) {
        // t_id_t_i; index of treated unit to match among treated units
        // t_id_i: index of treated unit to match among all units
        counter++;
        if (counter == 200) {
          counter = 0;
          Rcpp::checkUserInterrupt();
        }

        if (max(as<IntegerVector>(n_eligible[g_c])) == 0) {
          break;
        }

        t_id_i = ind_focal[t_id_t_i];

        p.increment();

        if (!eligible[t_id_i]) {
          continue;
        }

        k_total = 0;

        for (int gi : g_c) {
          k = find_control_mat(t_id_i,
                                treat,
                                ind_non_focal,
                                distance_mat.row(t_id_t_i),
                                eligible,
                                gi,
                                r,
                                mm.row(t_id_t_i),
                                ncc,
                                caliper_covs_mat,
                                caliper_covs,
                                caliper_dist,
                                use_exact,
                                exact,
                                aenc,
                                antiexact_covs);

          if (k.empty()) {
            if (r == 1) {
              k_total = 0;
              break;
            }
            continue;
          }

          matches_i[k_total] = k[0];
          k_total++;
        }

        if (k_total == 0) {
          eligible[t_id_i] = false;
          n_eligible[focal]--;
          continue;
        }

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
    int t_id_t_i;
    for (int t_id_t_i_ : ord) {
      // t_id_t_i; index of treated unit to match among treated units
      // t_id_i: index of treated unit to match among all units
      counter++;
      if (counter == 200) {
        counter = 0;
        Rcpp::checkUserInterrupt();
      }

      t_id_t_i = t_id_t_i_ - 1;

      t_id_i = ind_focal[t_id_t_i];

      p.increment();

      if (!eligible[t_id_i]) {
        continue;
      }

      k_total = 0;

      for (int gi : g_c) {
        k = find_control_mat(t_id_i,
                              treat,
                              ind_non_focal,
                              distance_mat.row(t_id_t_i),
                              eligible,
                              gi,
                              1,
                              mm.row(t_id_t_i),
                              ncc,
                              caliper_covs_mat,
                              caliper_covs,
                              caliper_dist,
                              use_exact,
                              exact,
                              aenc,
                              antiexact_covs,
                              ratio[t_id_t_i]);

        if (k.empty()) {
          k_total = 0;
          break;
        }

        for (int cc : k) {
          matches_i[k_total] = cc;
          k_total++;
        }
      }

      if (k_total == 0) {
        continue;
      }

      for (c = 0; c < k_total; c++) {
        mm(t_id_t_i, sum(!is_na(mm(t_id_t_i, _)))) = matches_i[c];
      }
    }
  }

  p.update(prog_length);

  mm = mm + 1;
  rownames(mm) = as<CharacterVector>(lab[ind_focal]);

  return mm;
}