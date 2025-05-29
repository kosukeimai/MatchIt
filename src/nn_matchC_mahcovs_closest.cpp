// [[Rcpp::depends(RcppProgress)]]
#include "eta_progress_bar.h"
#include "internal.h"
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
      ind_match[indt_tmp[i]] = i;
    }
  }

  IntegerVector ind_focal = indt[Range(indt_sep[focal], indt_sep[focal + 1] - 1)];

  std::vector<int> times_matched(n, 0);

  std::vector<int> times_matched_allowed(n, reuse_max);
  for (i = 0; i < nf; i++) {
    times_matched_allowed[ind_focal[i]] = ratio[i];
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

  //distance & caliper_dist
  bool use_caliper_dist = false;
  double caliper_dist;
  NumericVector distance;
  if (caliper_dist_.isNotNull() && distance_.isNotNull()) {
    distance = as<NumericVector>(distance_);
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

  // Matching variable: if any mah_covs equal to caliper_covs, use
  // that caliper_covs and caliper as matching variable
  NumericVector match_var;
  double match_var_caliper = R_PosInf;
  int n_mah_covs = mah_covs.ncol();
  if (ncc > 0) {
    double a;
    for (int mci = 0; match_var.size() == 0 && mci < n_mah_covs; mci++) {
      for (int cci = 0; match_var.size() == 0 && cci < ncc; cci++) {
        if (caliper_covs[cci] < 0) {
          continue;
        }

        a = get_affine_transformation(caliper_covs_mat.column(cci),
                                      mah_covs.column(mci));

        if (std::abs(a) > 1e-10) {
          caliper_covs_mat.column(cci) = mah_covs.column(mci);
          caliper_covs[cci] *= a;

          match_var = mah_covs.column(mci);
          match_var_caliper = caliper_covs[cci];
        }
      }
    }
  }

  if (match_var.size() == 0) {
    match_var = mah_covs.column(0);
  }

  IntegerVector ind_d_ord = o(match_var);
  ind_d_ord = ind_d_ord - 1; //location of each unit after sorting

  IntegerVector match_d_ord = o(ind_d_ord);
  match_d_ord = match_d_ord - 1;

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

  IntegerVector heap_ord;
  std::vector<int> k;
  k.reserve(1);
  R_xlen_t hi;

  std::function<bool(int, int)> cmp;
  if (close) {
    cmp = [&dist](const int& a, const int& b) {return dist[a] < dist[b];};
  }
  else {
    cmp = [&dist](const int& a, const int& b) {return dist[a] >= dist[b];};
  }

  IntegerVector::iterator ci;

  for (r = 1; r <= max_ratio; r++) {
    //Find closest control unit to each treated unit
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

      k = find_control_mahcovs(ti,
                               ind_d_ord,
                               match_d_ord,
                               match_var,
                               match_var_caliper,
                               treat,
                               distance,
                               eligible,
                               gi,
                               r,
                               mm.row(t_id_t_i),
                               mah_covs,
                               ncc,
                               caliper_covs_mat,
                               caliper_covs,
                               use_caliper_dist,
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
      dist.push_back(sum(pow(mah_covs.row(ti) - mah_covs.row(k[0]), 2.0)));
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

        k = find_control_mahcovs(t_id_i,
                                 ind_d_ord,
                                 match_d_ord,
                                 match_var,
                                 match_var_caliper,
                                 treat,
                                 distance,
                                 eligible,
                                 gi,
                                 r,
                                 mm.row(t_id_t_i),
                                 mah_covs,
                                 ncc,
                                 caliper_covs_mat,
                                 caliper_covs,
                                 use_caliper_dist,
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
        dist[hi] = sum(pow(mah_covs.row(t_id_i) - mah_covs.row(k[0]), 2.0));

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