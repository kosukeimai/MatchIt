// [[Rcpp::depends(RcppProgress)]]
#include "eta_progress_bar.h"
#include "internal.h"
using namespace Rcpp;

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

  //Next column to fill in each row of `mm`. Tracked rather than recomputed with
  //`sum(!is_na(mm(row, _)))`, which allocates twice for every match written.
  std::vector<int> mm_filled(mm.nrow(), 0);

  const CharacterVector lab = treat.names();

  //`base::order()`'s radix sort beats every C++ alternative measured here by 3-8x at
  //these sizes; see _dev/cpp-cleanup-notes.md. Looked up in the base environment
  //because `Function("order")` searches from the global environment, where a user
  //object of that name would mask it.
  Function o = Environment::base_env()["order"];

  //`as<>()` on a `Nullable` wraps the caller's SEXP rather than copying it, so every
  //object taken from an argument below is `const`. Writing through one of them would
  //modify the R object the caller passed in, and the change would outlive the call.

  //exact
  const bool use_exact = exact_.isNotNull();
  const IntegerVector exact = use_exact ? as<IntegerVector>(exact_) : IntegerVector(0);

  //distance & caliper_dist
  const bool use_caliper_dist = caliper_dist_.isNotNull() && distance_.isNotNull();
  const NumericVector distance = use_caliper_dist ? as<NumericVector>(distance_) : NumericVector(0);
  const double caliper_dist = use_caliper_dist ? as<double>(caliper_dist_) : R_PosInf;

  //caliper_covs
  const NumericVector caliper_covs = caliper_covs_.isNotNull() ? as<NumericVector>(caliper_covs_) : NumericVector(0);
  const NumericMatrix caliper_covs_mat = caliper_covs_.isNotNull() ? as<NumericMatrix>(caliper_covs_mat_) : NumericMatrix(0, 0);
  const int ncc = caliper_covs_mat.ncol();

  //antiexact
  const IntegerMatrix antiexact_covs = antiexact_covs_.isNotNull() ? as<IntegerMatrix>(antiexact_covs_) : IntegerMatrix(0, 0);
  const int aenc = antiexact_covs.ncol();

  //unit_id
  const bool use_unit_id = unit_id_.isNotNull();
  const IntegerVector unit_id = use_unit_id ? as<IntegerVector>(unit_id_) : IntegerVector(0);

  //Matching variable: when a caliper covariate is an affine transformation of one of
  //the `mah_covs` columns, sorting on that column lets the scan in
  //`find_control_mahcovs()` stop as soon as the caliper is exceeded. Only the caliper
  //is converted to that column's scale, and only for this function's own use; the
  //caliper is still enforced on its own scale by `caliper_covs_okay()`, and
  //`caliper_covs` and `caliper_covs_mat` belong to the caller and are left alone.
  const int n_mah_covs = mah_covs.ncol();
  int match_var_col = 0;
  double match_var_caliper = R_PosInf;
  bool match_var_found = false;

  for (int mci = 0; !match_var_found && mci < n_mah_covs; mci++) {
    for (int cci = 0; cci < ncc; cci++) {
      if (caliper_covs[cci] < 0) {
        continue;
      }

      double a = get_affine_transformation(caliper_covs_mat.column(cci),
                                           mah_covs.column(mci));

      if (std::abs(a) <= 1e-10) {
        continue;
      }

      match_var_col = mci;
      match_var_caliper = std::abs(a) * caliper_covs[cci];
      match_var_found = true;
      break;
    }
  }

  const NumericVector match_var = mah_covs.column(match_var_col);

  IntegerVector ind_d_ord = o(match_var);
  ind_d_ord = ind_d_ord - 1; //location of each unit after sorting

  //`ind_d_ord` is a permutation, so its order is just its inverse; computing that
  //directly avoids a second call into R
  IntegerVector match_d_ord(n);
  for (i = 0; i < n; i++) {
    match_d_ord[ind_d_ord[i]] = static_cast<int>(i);
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

  IntegerVector heap_ord;
  std::vector<int> k;
  k.reserve(1);
  R_xlen_t hi;

  //One lambda rather than two wrapped in a `std::function`, so the comparison can
  //be inlined into `std::lower_bound()` below
  auto cmp = [&dist, close](const int& a, const int& b) {
    if (close) {
      return dist[a] < dist[b];
    }

    return dist[a] >= dist[b];
  };

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
      dist.push_back(euc_dist_sq(mah_covs, ti, k[0]));
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
        dist[hi] = euc_dist_sq(mah_covs, t_id_i, k[0]);

        // Find new position of pair in heap
        ci = std::lower_bound(heap_ord.begin() + i, heap_ord.end(), hi, cmp);

        if (ci != heap_ord.begin() + i) {
          std::rotate(heap_ord.begin() + i, heap_ord.begin() + i + 1, ci);
        }

        continue;
      }

      mm(t_id_t_i, mm_filled[t_id_t_i]++) = c_id_i;

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
