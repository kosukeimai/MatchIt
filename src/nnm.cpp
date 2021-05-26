// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
IntegerMatrix nn_matchC(const IntegerVector& treat_,
                        const IntegerVector& ord_,
                        const IntegerVector& ratio,
                        const int& max_rat,
                        const bool& replace,
                        const LogicalVector& discarded,
                        const Nullable<NumericVector>& distance_ = R_NilValue,
                        const Nullable<NumericMatrix>& distance_mat_ = R_NilValue,
                        const Nullable<IntegerVector>& exact_ = R_NilValue,
                        const Nullable<double>& caliper_dist_ = R_NilValue,
                        const Nullable<NumericVector>& caliper_covs_ = R_NilValue,
                        const Nullable<NumericMatrix>& calcovs_covs_mat_ = R_NilValue,
                        const Nullable<NumericMatrix>& mah_covs_ = R_NilValue,
                        const Nullable<NumericMatrix>& mahSigma_inv_ = R_NilValue,
                        const Nullable<IntegerMatrix>& antiexact_covs_ = R_NilValue,
                        const bool& disl_prog = false)
  {

  // Initialize

  NumericVector distance, caliper_covs;
  double caliper_dist;
  NumericMatrix distance_mat, calcovs_covs_mat, mah_covs, mahSigma_inv, mah_covs_c;
  IntegerMatrix antiexact_covs;
  IntegerVector exact, exact_c, antiexact_col;

  Environment pkg = Environment::namespace_env("stats");
  Function mah = pkg["mahalanobis"];

  bool use_dist_mat = false;
  bool use_exact = false;
  bool use_caliper_dist = false;
  bool use_caliper_covs = false;
  bool use_mah_covs = false;
  bool use_antiexact = false;

  // Info about original treat
  int n_ = treat_.size();
  IntegerVector ind_ = Range(0, n_ - 1);
  IntegerVector ind1_ = ind_[treat_ == 1];
  IntegerVector ind0_ = ind_[treat_ == 0];
  int n1_ = ind1_.size();
  int n0_ = n_ - n1_;

  // Output matrix with sample indices of C units
  IntegerMatrix mm(n1_, max_rat);
  mm.fill(NA_INTEGER);
  rownames(mm) = as<CharacterVector>(treat_.names())[ind1_];

  // Store who has been matched
  LogicalVector matched = clone(discarded);

  // After discarding

  IntegerVector ind = ind_[!discarded];
  IntegerVector treat = treat_[!discarded];
  IntegerVector ind0 = ind[treat == 0];
  int n0 = ind0.size();

  int t, t_ind, min_ind, c_chosen, num_eligible, cal_len, t_rat, n_anti;
  double dt, cal_var_t;

  NumericVector cal_var, cal_diff, ps_diff, diff, dist_t, mah_covs_t, mah_covs_col,
                match_distance(n0);

  IntegerVector c_eligible(n0), indices(n0);
  LogicalVector finite_match_distance(n0);

  if (distance_.isNotNull()) {
    distance = distance_;
  }
  if (exact_.isNotNull()) {
    exact = exact_;
    use_exact = true;
  }
  if (caliper_dist_.isNotNull()) {
    caliper_dist = as<double>(caliper_dist_);
    use_caliper_dist = true;
    ps_diff = NumericVector(n_);
  }
  if (caliper_covs_.isNotNull()) {
    caliper_covs = caliper_covs_;
    use_caliper_covs = true;
    cal_len = caliper_covs.size();
    cal_diff = NumericVector(n0);
  }
  if (calcovs_covs_mat_.isNotNull()) {
    calcovs_covs_mat = as<NumericMatrix>(calcovs_covs_mat_);
  }
  if (mah_covs_.isNotNull()) {
    mah_covs = as<NumericMatrix>(mah_covs_);
    use_mah_covs = true;
  } else {
    if (distance_mat_.isNotNull()) {
      distance_mat = as<NumericMatrix>(distance_mat_);

      // IntegerVector ind0_ = ind_[treat_ == 0];
      NumericVector dist_t(n0_);
      use_dist_mat = true;
    }
    ps_diff = NumericVector(n_);
  }
  if (mahSigma_inv_.isNotNull()) {
    mahSigma_inv = as<NumericMatrix>(mahSigma_inv_);
  }
  if (antiexact_covs_.isNotNull()) {
    antiexact_covs = as<IntegerMatrix>(antiexact_covs_);
    n_anti = antiexact_covs.ncol();
    use_antiexact = true;
  }

  bool ps_diff_assigned;

  //progress bar
  int prog_length;
  if (replace) prog_length = n1_ + 1;
  else prog_length = max_rat*n1_ + 1;
  Progress p(prog_length, disl_prog);

  //Counters
  int rat, i, x, j, a, k;
  k = -1;

  //Matching
  for (rat = 0; rat < max_rat; ++rat) {
    for (i = 0; i < n1_; ++i) {

      k++;
      if (k % 500 == 0) Rcpp::checkUserInterrupt();

      p.increment();

      if (all(as<LogicalVector>(matched[ind0])).is_true()) {
        break;
      }

      t = ord_[i] - 1;   // index among treated
      t_ind = ind1_[t]; // index among sample

      if (matched[t_ind]) {
        continue;
      }

      //Check if unit has enough matches
      t_rat = ratio[t];

      if (t_rat < rat + 1) {
        continue;
      }

      c_eligible = ind0; // index among sample

      c_eligible = c_eligible[!as<LogicalVector>(matched[c_eligible])];

      if (use_exact) {
        exact_c = exact[c_eligible];
        c_eligible = c_eligible[exact_c == exact[t_ind]];
      }

      if (c_eligible.size() == 0) {
        continue;
      }

      if (use_antiexact) {
        for (a = 0; a < n_anti; ++a) {
          antiexact_col = antiexact_covs(_, a);
          antiexact_col = antiexact_col[c_eligible];
          c_eligible = antiexact_col[antiexact_col != antiexact_col[t_ind]];
        }
      }

      if (c_eligible.size() == 0) {
        continue;
      }

      ps_diff_assigned = false;

      if (use_caliper_dist) {
        if (use_dist_mat) {
          dist_t = distance_mat.row(t);
          diff = dist_t[match(c_eligible, ind0_) - 1];
        } else {
          dt = distance[t_ind];
          diff = Rcpp::abs(as<NumericVector>(distance[c_eligible]) - dt);
        }

        ps_diff[c_eligible] = diff;
        ps_diff_assigned = true;

        c_eligible = c_eligible[diff <= caliper_dist];

        if (c_eligible.size() == 0) {
          continue;
        }
      }

      if (use_caliper_covs) {
        for (x = 0; (x < cal_len) && c_eligible.size() > 0; ++x) {
          cal_var = calcovs_covs_mat( _ , x );

          cal_var_t = cal_var[t_ind];

          diff = Rcpp::abs(as<NumericVector>(cal_var[c_eligible]) - cal_var_t);

          cal_diff = diff;

          c_eligible = c_eligible[cal_diff <= caliper_covs[x]];
        }

        if (c_eligible.size() == 0) {
          continue;
        }
      }

      //Compute distances among eligible
      num_eligible = c_eligible.size();

      //If replace and few eligible controls, assign all and move on
      if (replace && (num_eligible <= t_rat)) {
        for (j = 0; j < num_eligible; ++j) {
          mm( t , j ) = c_eligible[j] + 1;
        }
        continue;
      }

      if (use_mah_covs) {
        mah_covs_c = NumericMatrix(num_eligible, mah_covs.ncol());
        for (j = 0; j < mah_covs.ncol(); ++j) {
          mah_covs_col = mah_covs.column(j);
          mah_covs_c(_,j) = as<NumericVector>(mah_covs_col[c_eligible]);
        }
        mah_covs_t = mah_covs( t_ind , _ );
        match_distance = sqrt(as<NumericVector>(mah(mah_covs_c, mah_covs_t, mahSigma_inv, true))); //mahalanobis in R

      } else if (ps_diff_assigned) {
        match_distance = ps_diff[c_eligible]; //c_eligible might have shrunk since previous assignment
      } else if (use_dist_mat) {
        dist_t = distance_mat.row(t);
        match_distance = dist_t[match(c_eligible, ind0_) - 1];
      } else {
        dt = distance[t_ind];
        match_distance = Rcpp::abs(as<NumericVector>(distance[c_eligible]) - dt);
      }

      //Remove infinite distances
      finite_match_distance = is_finite(match_distance);
      c_eligible = c_eligible[finite_match_distance];
      if (c_eligible.size() == 0) {
        continue;
      }
      match_distance = match_distance[finite_match_distance];

      if (replace) {
        //When matching w/ replacement, get t_rat closest control units
        indices = Range(0, num_eligible - 1);

        std::partial_sort(indices.begin(), indices.begin() + t_rat, indices.end(),
                          [&match_distance](int k, int j) {return match_distance[k] < match_distance[j];});

        for (j = 0; j < t_rat; ++j) {
          min_ind = indices[j];
          mm( t , j ) = c_eligible[min_ind] + 1;
        }
      }
      else {
        min_ind = which_min(match_distance);
        c_chosen = c_eligible[min_ind];

        mm( t , rat ) = c_chosen + 1; // + 1 because C indexing starts at 0 but mm is sent to R

        matched[c_chosen] = true;
      }
    }

    if (replace) break;
  }

  p.update(prog_length);

  return mm;
}
