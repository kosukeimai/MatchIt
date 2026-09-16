// [[Rcpp::depends(RcppProgress)]]
#include "eta_progress_bar.h"
#include "internal.h"
using namespace Rcpp;

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

  IntegerVector ind_d_ord = o(distance);
  ind_d_ord = ind_d_ord - 1; //location of each unit after sorting

  //`ind_d_ord` is a permutation, so its order is just its inverse; computing that
  //directly avoids a second call into R
  IntegerVector match_d_ord(n);
  for (i = 0; i < n; i++) {
    match_d_ord[ind_d_ord[i]] = static_cast<int>(i);
  }

  IntegerVector last_control(g);
  last_control.fill(n - 1);
  IntegerVector first_control(g);
  first_control.fill(0);

  //`as<>()` on a `Nullable` wraps the caller's SEXP rather than copying it, so every
  //object taken from an argument below is `const`. Writing through one of them would
  //modify the R object the caller passed in, and the change would outlive the call.

  //exact
  const bool use_exact = exact_.isNotNull();
  const IntegerVector exact = use_exact ? as<IntegerVector>(exact_) : IntegerVector(0);

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

  //caliper_dist. Not `const`: the loop below may tighten it.
  double caliper_dist = caliper_dist_.isNotNull() ? as<double>(caliper_dist_) : max_finite(distance) - min_finite(distance) + 1;

  //A caliper on a covariate that is an affine transformation of `distance` is
  //equivalent to a caliper on `distance` itself, which the sorted scan can stop early
  //on; the covariate caliper stays in force either way.
  for (int cci = 0; cci < ncc; cci++) {
    double a = get_affine_transformation(caliper_covs_mat.column(cci),
                                         distance);

    if (std::abs(a) <= 1e-10) {
      continue;
    }

    //`std::abs()` because a negative `a` would otherwise flip the sign of the
    //caliper, and a negative caliper means the opposite of a positive one
    double caliper_dist_cci = std::abs(a) * caliper_covs[cci];

    if (caliper_dist_.isNull() ||
        (caliper_covs[cci] >= 0 && caliper_dist > caliper_dist_cci) ||
        (caliper_covs[cci] < 0 && caliper_dist < caliper_dist_cci)) {
      caliper_dist = caliper_dist_cci;
    }
  }

  //storing closeness
  std::vector<int> t_id, c_id;
  std::vector<double> dist;
  t_id.reserve(n_eligible[focal]);
  c_id.reserve(n_eligible[focal]);
  dist.reserve(n_eligible[focal]);

  gi = 0;

  update_first_and_last_control(first_control,
                                last_control,
                                ind_d_ord,
                                eligible,
                                treat,
                                gi);

  IntegerVector ck_;

  int c_id_i, t_id_t_i, t_id_i;

  int counter = 0;
  int r = 1;

  IntegerVector heap_ord;
  std::vector<int> k(1);
  R_xlen_t hi;

  //progress bar
  R_xlen_t prog_length = sum(ratio) + 1;
  ETAProgressBar pb;
  Progress p(prog_length, disl_prog, pb);

  IntegerVector::iterator ci;

  //One lambda rather than two wrapped in a `std::function`, so the comparison can
  //be inlined into `std::lower_bound()` below
  auto cmp = [&dist, close](const int& a, const int& b) {
    if (close) {
      return dist[a] < dist[b];
    }

    return dist[a] >= dist[b];
  };

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

      k = find_control_vec(ti,
                           ind_d_ord,
                           match_d_ord,
                           treat,
                           distance,
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
                           antiexact_covs,
                           first_control,
                           last_control);

      if (k.empty()) {
        eligible[ti] = false;
        n_eligible[focal]--;
        continue;
      }

      t_id.push_back(ti);
      c_id.push_back(k[0]);
      dist.push_back(std::abs(distance[ti] - distance[k[0]]));
    }

    nf = dist.size();

    //Order the list
    heap_ord = o(dist, _["decreasing"] = !close);
    heap_ord = heap_ord - 1;

    i = 0;

    //Go through ordered list and assign matches, re-matching when necessary
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
        update_first_and_last_control(first_control,
                                      last_control,
                                      ind_d_ord,
                                      eligible,
                                      treat,
                                      gi);

        k = find_control_vec(t_id_i,
                             ind_d_ord,
                             match_d_ord,
                             treat,
                             distance,
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
                             antiexact_covs,
                             first_control,
                             last_control,
                             1,
                             c_id_i);

        //If no matches...
        if (k.empty()) {
          eligible[t_id_i] = false;
          n_eligible[focal]--;
          continue;
        }

        c_id[hi] = k[0];
        dist[hi] = std::abs(distance[t_id_i] - distance[k[0]]);

        // Find new position of pair in heap
        ci = std::lower_bound(heap_ord.begin() + i, heap_ord.end(), hi,
                              cmp);

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
