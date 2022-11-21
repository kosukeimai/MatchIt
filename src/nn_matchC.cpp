// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
IntegerMatrix nn_matchC(const IntegerVector& treat_,
                        const IntegerVector& ord_,
                        const IntegerVector& ratio,
                        const LogicalVector& discarded,
                        const int& reuse_max,
                        const Nullable<NumericVector>& distance_ = R_NilValue,
                        const Nullable<NumericMatrix>& distance_mat_ = R_NilValue,
                        const Nullable<IntegerVector>& exact_ = R_NilValue,
                        const Nullable<double>& caliper_dist_ = R_NilValue,
                        const Nullable<NumericVector>& caliper_covs_ = R_NilValue,
                        const Nullable<NumericMatrix>& caliper_covs_mat_ = R_NilValue,
                        const Nullable<NumericMatrix>& mah_covs_ = R_NilValue,
                        const Nullable<IntegerMatrix>& antiexact_covs_ = R_NilValue,
                        const Nullable<IntegerVector>& unit_id_ = R_NilValue,
                        const bool& disl_prog = false)
  {

  // Initialize

  NumericVector distance, caliper_covs;
  double caliper_dist;
  NumericMatrix distance_mat, caliper_covs_mat, mah_covs, mah_covs_c;
  IntegerMatrix antiexact_covs;
  IntegerVector exact, exact_c, antiexact_col, unit_id, units_with_id_of_chosen_unit;
  int id_of_chosen_unit;

  bool use_dist_mat = false;
  bool use_exact = false;
  bool use_caliper_dist = false;
  bool use_caliper_covs = false;
  bool use_mah_covs = false;
  bool use_antiexact = false;
  bool use_reuse_max = false;

  // Info about original treat
  int n_ = treat_.size();
  IntegerVector ind_ = Range(0, n_ - 1);
  IntegerVector ind1_ = ind_[treat_ == 1];
  IntegerVector ind0_ = ind_[treat_ == 0];
  int n1_ = ind1_.size();
  int n0_ = n_ - n1_;
  CharacterVector lab = treat_.names();

  // Output matrix with sample indices of C units
  int max_rat = max(ratio);
  IntegerMatrix mm(n1_, max_rat);
  mm.fill(NA_INTEGER);
  rownames(mm) = lab[ind1_];

  // Store who has been matched
  IntegerVector matched = rep(0, n_);
  matched[discarded] = n1_; //discarded are unmatchable

  // After discarding

  IntegerVector ind = ind_[!discarded];
  IntegerVector treat = treat_[!discarded];
  IntegerVector ind0 = ind[treat == 0];
  int n0 = ind0.size();

  int t, t_ind, min_ind, c_chosen, num_eligible, cal_len, t_rat, n_anti, antiexact_t;
  double dt, cal_var_t;

  NumericVector cal_var, cal_diff, ps_diff, diff, dist_t, mah_covs_t, mah_covs_row,
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
  if (caliper_covs_mat_.isNotNull()) {
    caliper_covs_mat = as<NumericMatrix>(caliper_covs_mat_);
  }
  if (mah_covs_.isNotNull()) {
    mah_covs = as<NumericMatrix>(mah_covs_);
    NumericVector mah_covs_row(mah_covs.ncol());
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
  if (antiexact_covs_.isNotNull()) {
    antiexact_covs = as<IntegerMatrix>(antiexact_covs_);
    n_anti = antiexact_covs.ncol();
    use_antiexact = true;
  }
  if (reuse_max < n1_) {
    use_reuse_max = true;
  }
  if (unit_id_.isNotNull()) {
    unit_id = as<IntegerVector>(unit_id_);
  }
  else {
    unit_id = ind_;
  }

  bool ps_diff_assigned;

  //progress bar
  int prog_length;
  if (!use_reuse_max) prog_length = n1_ + 1;
  else prog_length = max_rat*n1_ + 1;
  Progress p(prog_length, disl_prog);

  //Counters
  int rat, i, x, j, j_, a, k;
  k = -1;

  //Matching
  for (rat = 0; rat < max_rat; ++rat) {
    for (i = 0; i < n1_; ++i) {

      k++;
      if (k % 500 == 0) Rcpp::checkUserInterrupt();

      p.increment();

      if (all(as<IntegerVector>(matched[ind0]) >= reuse_max).is_true()){
        break;
      }

      t = ord_[i] - 1;   // index among treated
      t_ind = ind1_[t];  // index among sample

      // Skip discarded units (only discarded treated units have matched > 0)
      if (matched[t_ind] > 0) {
        continue;
      }

      //Check if unit has enough matches
      t_rat = ratio[t];

      if (t_rat < rat + 1) {
        continue;
      }

      c_eligible = ind0; // index among sample

      //Make ineligible control units that have been matched too many times
      c_eligible = c_eligible[as<IntegerVector>(matched[c_eligible]) < reuse_max];

      //Prevent control units being matched to same treated unit again
      if (rat > 0) {
        c_eligible = as<IntegerVector>(c_eligible[!in(c_eligible, na_omit(mm.row(t)) - 1)]);
      }

      if (use_exact) {
        exact_c = exact[c_eligible];
        c_eligible = c_eligible[exact_c == exact[t_ind]];
      }

      if (c_eligible.size() == 0) {
        continue;
      }

      if (use_antiexact) {
        for (a = 0; a < n_anti && c_eligible.size() > 0; ++a) {
          antiexact_col = antiexact_covs(_, a);
          antiexact_t = antiexact_col[t_ind];
          antiexact_col = antiexact_col[c_eligible];
          c_eligible = c_eligible[antiexact_col != antiexact_t];
        }
        if (c_eligible.size() == 0) {
          continue;
        }
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
          cal_var = caliper_covs_mat( _ , x );

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
      if (!use_reuse_max && (num_eligible <= t_rat)) {
        for (j = 0; j < num_eligible; ++j) {
          mm( t , j ) = c_eligible[j] + 1;
        }
        continue;
      }

      if (use_mah_covs) {

        match_distance = rep(0.0, num_eligible);
        mah_covs_t = mah_covs.row(t_ind);

        for (j = 0; j < num_eligible; j++) {
          j_ = c_eligible[j];
          mah_covs_row = mah_covs.row(j_);
          match_distance[j] = sqrt(sum(pow(mah_covs_t - mah_covs_row, 2.0)));
        }

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
      num_eligible = c_eligible.size();
      if (num_eligible == 0) {
        continue;
      }
      match_distance = match_distance[finite_match_distance];

      if (!use_reuse_max) {
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

        id_of_chosen_unit = unit_id[c_chosen];
        units_with_id_of_chosen_unit = ind_[unit_id == id_of_chosen_unit];
        for (j = 0; j < units_with_id_of_chosen_unit.size(); ++j) {
          matched[units_with_id_of_chosen_unit[j]] = matched[units_with_id_of_chosen_unit[j]] + 1;
        }
      }
    }

    if (!use_reuse_max) break;
  }

  p.update(prog_length);

  return mm;
}
