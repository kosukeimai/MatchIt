#ifndef INTERNAL_H
#define INTERNAL_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>
#include <vector>

//No `using namespace Rcpp;` here: this header is included by every translation unit
//in the package, and each of them has its own `using` directive.

Rcpp::IntegerVector tabulateC_(const Rcpp::IntegerVector& bins,
                               int nbins = 0);

Rcpp::IntegerVector which(const Rcpp::LogicalVector& x);

int recode_focal(int focal,
                 const Rcpp::IntegerVector& unique_treat);

std::vector<int> find_control_vec(int t_id,
                                  const Rcpp::IntegerVector& ind_d_ord,
                                  const Rcpp::IntegerVector& match_d_ord,
                                  const Rcpp::IntegerVector& treat,
                                  const Rcpp::NumericVector& distance,
                                  const Rcpp::LogicalVector& eligible,
                                  int gi,
                                  int r,
                                  const Rcpp::IntegerVector& mm_rowi_,
                                  int ncc,
                                  const Rcpp::NumericMatrix& caliper_covs_mat,
                                  const Rcpp::NumericVector& caliper_covs,
                                  double caliper_dist,
                                  bool use_exact,
                                  const Rcpp::IntegerVector& exact,
                                  int aenc,
                                  const Rcpp::IntegerMatrix& antiexact_covs,
                                  const Rcpp::IntegerVector& first_control,
                                  const Rcpp::IntegerVector& last_control,
                                  int ratio = 1,
                                  int prev_start = -1);

std::vector<int> find_control_mahcovs(int t_id,
                                      const Rcpp::IntegerVector& ind_d_ord,
                                      const Rcpp::IntegerVector& match_d_ord,
                                      const Rcpp::NumericVector& match_var,
                                      double match_var_caliper,
                                      const Rcpp::IntegerVector& treat,
                                      const Rcpp::NumericVector& distance,
                                      const Rcpp::LogicalVector& eligible,
                                      int gi,
                                      int r,
                                      const Rcpp::IntegerVector& mm_rowi,
                                      const Rcpp::NumericMatrix& mah_covs,
                                      int ncc,
                                      const Rcpp::NumericMatrix& caliper_covs_mat,
                                      const Rcpp::NumericVector& caliper_covs,
                                      bool use_caliper_dist,
                                      double caliper_dist,
                                      bool use_exact,
                                      const Rcpp::IntegerVector& exact,
                                      int aenc,
                                      const Rcpp::IntegerMatrix& antiexact_covs,
                                      int ratio = 1);

std::vector<int> find_control_mat(int t_id,
                                  const Rcpp::IntegerVector& treat,
                                  const Rcpp::IntegerVector& ind_non_focal,
                                  const Rcpp::NumericVector& distance_mat_row_i,
                                  const Rcpp::LogicalVector& eligible,
                                  int gi,
                                  int r,
                                  const Rcpp::IntegerVector& mm_rowi,
                                  int ncc,
                                  const Rcpp::NumericMatrix& caliper_covs_mat,
                                  const Rcpp::NumericVector& caliper_covs,
                                  double caliper_dist,
                                  bool use_exact,
                                  const Rcpp::IntegerVector& exact,
                                  int aenc,
                                  const Rcpp::IntegerMatrix& antiexact_covs,
                                  int ratio = 1);

double euc_dist_sq(const Rcpp::NumericMatrix& x,
                   int i,
                   int j);

//`ids` is taken by value so the early returns can move it rather than copy
std::vector<int> take_closest(std::vector<int> ids,
                              const std::vector<double>& dists,
                              int ratio);

bool antiexact_okay(int aenc,
                    int i,
                    int j,
                    const Rcpp::IntegerMatrix& antiexact_covs);

bool caliper_covs_okay(int ncc,
                       int i,
                       int j,
                       const Rcpp::NumericMatrix& caliper_covs_mat,
                       const Rcpp::NumericVector& caliper_covs);

bool caliper_dist_okay(bool use_caliper_dist,
                       int i,
                       int j,
                       const Rcpp::NumericVector& distance,
                       double caliper_dist);

bool mm_okay(int r,
             int i,
             const Rcpp::IntegerVector& mm_rowi);

bool exact_okay(bool use_exact,
                int i,
                int j,
                const Rcpp::IntegerVector& exact);

double max_finite(const Rcpp::NumericVector& x);

double min_finite(const Rcpp::NumericVector& x);

//`first_control` and `last_control` are taken by value on purpose: the copies share
//the caller's SEXP, which is how the updates below reach the caller.
void update_first_and_last_control(Rcpp::IntegerVector first_control,
                                   Rcpp::IntegerVector last_control,
                                   const Rcpp::IntegerVector& ind_d_ord,
                                   const Rcpp::LogicalVector& eligible,
                                   const Rcpp::IntegerVector& treat,
                                   int gi);

double get_affine_transformation(const Rcpp::NumericVector& x,
                                 const Rcpp::NumericVector& y,
                                 double tol = 1e-9);

#endif
