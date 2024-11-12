#ifndef INTERNAL_H
#define INTERNAL_H

#include <Rcpp.h>
#include <algorithm>
#include <functional>
#include <cmath>
#include <utility>
#include <tuple>
using namespace Rcpp;

IntegerVector tabulateC_(const IntegerVector& bins,
                         const int& nbins = 0);

IntegerVector which(const LogicalVector& x);

std::vector<int> find_control_vec(const int& t_id,
                                  const IntegerVector& ind_d_ord,
                                  const IntegerVector& match_d_ord,
                                  const IntegerVector& treat,
                                  const NumericVector& distance,
                                  const LogicalVector& eligible,
                                  const int& gi,
                                  const int& r,
                                  const IntegerVector& mm_rowi_,
                                  const int& ncc,
                                  const NumericMatrix& caliper_covs_mat,
                                  const NumericVector& caliper_covs,
                                  const double& caliper_dist,
                                  const bool& use_exact,
                                  const IntegerVector& exact,
                                  const int& aenc,
                                  const IntegerMatrix& antiexact_covs,
                                  const IntegerVector& first_control,
                                  const IntegerVector& last_control,
                                  const int& ratio = 1,
                                  const int& prev_start = -1);

std::vector<int> find_control_mahcovs(const int& t_id,
                                      const IntegerVector& ind_d_ord,
                                      const IntegerVector& match_d_ord,
                                      const NumericVector& match_var,
                                      const double& match_var_caliper,
                                      const IntegerVector& treat,
                                      const NumericVector& distance,
                                      const LogicalVector& eligible,
                                      const int& gi,
                                      const int& r,
                                      const IntegerVector& mm_rowi,
                                      const NumericMatrix& mah_covs,
                                      const int& ncc,
                                      const NumericMatrix& caliper_covs_mat,
                                      const NumericVector& caliper_covs,
                                      const bool& use_caliper_dist,
                                      const double& caliper_dist,
                                      const bool& use_exact,
                                      const IntegerVector& exact,
                                      const int& aenc,
                                      const IntegerMatrix& antiexact_covs,
                                      const int& ratio = 1);

std::vector<int> find_control_mat(const int& t_id,
                                  const IntegerVector& treat,
                                  const IntegerVector& ind_non_focal,
                                  const NumericVector& distance_mat_row_i,
                                  const LogicalVector& eligible,
                                  const int& gi,
                                  const int& r,
                                  const IntegerVector& mm_rowi,
                                  const int& ncc,
                                  const NumericMatrix& caliper_covs_mat,
                                  const NumericVector& caliper_covs,
                                  const double& caliper_dist,
                                  const bool& use_exact,
                                  const IntegerVector& exact,
                                  const int& aenc,
                                  const IntegerMatrix& antiexact_covs,
                                  const int& ratio = 1);

bool antiexact_okay(const int& aenc,
                    const int& i,
                    const int& j,
                    const IntegerMatrix& antiexact_covs);

bool caliper_covs_okay(const int& ncc,
                       const int& i,
                       const int& j,
                       const NumericMatrix& caliper_covs_mat,
                       const NumericVector& caliper_covs);

bool caliper_dist_okay(const bool& use_caliper_dist,
                       const int& i,
                       const int& j,
                       const NumericVector& distance,
                       const double& caliper_dist);

bool mm_okay(const int& r,
             const int& i,
             const IntegerVector& mm_rowi);

bool exact_okay(const bool& use_exact,
                const int& i,
                const int& j,
                const IntegerVector& exact);

double max_finite(const NumericVector& x);

double min_finite(const NumericVector& x);

void update_first_and_last_control(IntegerVector first_control,
                                   IntegerVector last_control,
                                   const IntegerVector& ind_d_ord,
                                   const LogicalVector& eligible,
                                   const IntegerVector& treat,
                                   const int& gi);

#endif
