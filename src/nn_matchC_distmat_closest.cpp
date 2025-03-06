// [[Rcpp::depends(RcppProgress)]]
#include "eta_progress_bar.h"
#include "internal.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
IntegerMatrix nn_matchC_distmat_closest(const IntegerVector& treat,
                                        const IntegerVector& ratio,
                                        const LogicalVector& discarded,
                                        const int& reuse_max,
                                        const NumericMatrix& distance_mat,
                                        const Nullable<IntegerMatrix>& exact_ = R_NilValue,
                                        const Nullable<double>& caliper_dist_ = R_NilValue,
                                        const Nullable<NumericVector>& caliper_covs_ = R_NilValue,
                                        const Nullable<NumericMatrix>& caliper_covs_mat_ = R_NilValue,
                                        const Nullable<IntegerMatrix>& antiexact_covs_ = R_NilValue,
                                        const Nullable<IntegerVector>& unit_id_ = R_NilValue,
                                        const bool& close = true,
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
    nt[gi] = std::count(treat.begin(), treat.end(), gi);
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

  IntegerVector ind_non_focal = which(treat != focal);
  IntegerVector ind_focal = which(treat == focal);

  for (i = 0; i < n - nf; i++) {
    ind_match[ind_non_focal[i]] = i;
  }

  for (i = 0; i < nf; i++) {
    ind_match[ind_focal[i]] = i;
  }

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

  Function o("order");

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

  //unit_id
  IntegerVector unit_id;
  bool use_unit_id = false;
  if (unit_id_.isNotNull()) {
    unit_id = as<IntegerVector>(unit_id_);
    use_unit_id = true;
  }

  //storing closeness
  std::vector<int> t_id, c_id;
  std::vector<double> dist;
  t_id.reserve(n_eligible[focal]);
  c_id.reserve(n_eligible[focal]);
  dist.reserve(n_eligible[focal]);

  //progress bar
  R_xlen_t prog_length = n_eligible[focal] + sum(ratio) + 1;
  ETAProgressBar pb;
  Progress p(prog_length, disl_prog, pb);

  gi = 0;

  IntegerVector ck_;

  int c_id_i, t_id_t_i, t_id_i;

  int counter = 0;
  int r = 1;

  IntegerVector heap_ord(n_eligible[focal]);
  std::vector<int> k;
  k.reserve(1);
  R_xlen_t hi;

  IntegerVector::iterator ci;

  std::function<bool(int, int)> cmp;
  if (close) {
    cmp = [&dist](const int& a, const int& b) {return dist[a] < dist[b];};
  }
  else {
    cmp = [&dist](const int& a, const int& b) {return dist[a] >= dist[b];};
  }

  for (r = 1; r <= max_ratio; r++) {
    for (int ti : ind_focal) {
      if (!eligible[ti]) {
        continue;
      }

      counter++;
      if (counter == 200) {
        counter = 0;
        Rcpp::checkUserInterrupt();
      }

      t_id_t_i = ind_match[ti];

      k = find_control_mat(ti,
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

      p.increment();

      if (k.empty()) {
        eligible[ti] = false;
        n_eligible[focal]--;
        continue;
      }

      t_id.push_back(ti);
      c_id.push_back(k[0]);
      dist.push_back(distance_mat(t_id_t_i, ind_match[k[0]]));
    }

    nf = dist.size();

    //Order the list
    heap_ord = o(dist, _["decreasing"] = !close);
    heap_ord = heap_ord - 1;

    i = 0;
    while (min(n_eligible) > 0 && i < nf) {
      counter++;
      if (counter == 200) {
        counter = 0;
        Rcpp::checkUserInterrupt();
      }

      hi = heap_ord[i];

      t_id_i = t_id[hi];

      if (!eligible[t_id_i]) {
        i++;
        continue;
      }

      t_id_t_i = ind_match[t_id_i];

      c_id_i = c_id[hi];

      if (!eligible[c_id_i]) {
        // If control isn't eligible, find new control and try again

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

        //If no matches...
        if (k.empty()) {
          eligible[t_id_i] = false;
          n_eligible[focal]--;
          continue;
        }

        c_id[hi] = k[0];
        dist[hi] = distance_mat(t_id_t_i, ind_match[k[0]]);

        // Find new position of pair in heap
        ci = std::lower_bound(heap_ord.begin() + i, heap_ord.end(), hi, cmp);

        if (ci != heap_ord.begin() + i) {
          std::rotate(heap_ord.begin() + i, heap_ord.begin() + i + 1, ci);
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

      i++;
    }

    t_id.clear();
    c_id.clear();
    dist.clear();
  }

  p.update(prog_length);

  mm = mm + 1;
  rownames(mm) = as<CharacterVector>(lab[ind_focal]);

  return mm;
}