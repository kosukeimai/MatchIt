#ifndef INTERNAL_H
#define INTERNAL_H

#include <Rcpp.h>
using namespace Rcpp;

IntegerVector tabulateC_(const IntegerVector& bins,
                         const int& nbins = 0);

IntegerVector which(const LogicalVector& x);

int find_both(const int& t_id,
              const IntegerVector& ind_d_ord,
              const IntegerVector& match_d_ord,
              const IntegerVector& treat,
              const NumericVector& distance,
              const LogicalVector& eligible,
              const int& gi,
              const int& r,
              const IntegerVector& mm_ordi,
              const int& ncc,
              const NumericMatrix& caliper_covs_mat,
              const NumericVector& caliper_covs,
              const double& caliper_dist,
              const bool& use_exact,
              const IntegerVector& exact,
              const int& aenc,
              const IntegerMatrix& antiexact_covs,
              const IntegerVector& first_control,
              const IntegerVector& last_control);

int find_lr(const int& prev_match,
            const int& t_id,
            const IntegerVector& ind_d_ord,
            const IntegerVector& match_d_ord,
            const IntegerVector& treat,
            const NumericVector& distance,
            const LogicalVector& eligible,
            const int& gi,
            const int& ncc,
            const NumericMatrix& caliper_covs_mat,
            const NumericVector& caliper_covs,
            const double& caliper_dist,
            const bool& use_exact,
            const IntegerVector& exact,
            const int& aenc,
            const IntegerMatrix& antiexact_covs,
            const IntegerVector& first_control,
            const IntegerVector& last_control);

bool antiexact_okay(const int& aenc,
                    const int& i,
                    const int& j,
                    const IntegerMatrix& antiexact_covs);

bool caliper_covs_okay(const int& ncc,
                       const int& i,
                       const int& j,
                       const NumericMatrix& caliper_covs_mat,
                       const NumericVector& caliper_covs);

bool mm_okay(const int& r,
             const int& i,
             const IntegerVector& mm_rowi);

bool exact_okay(const bool& use_exact,
                const int& i,
                const int& j,
                const IntegerVector& exact);

void swap_pos(IntegerVector x,
              const int& a,
              const int& b);

double max_finite(const NumericVector& x);

double min_finite(const NumericVector& x);

void update_first_and_last_control(IntegerVector first_control,
                                   IntegerVector last_control,
                                   const IntegerVector& ind_d_ord,
                                   const LogicalVector& eligible,
                                   const IntegerVector& treat,
                                   const int& gi);

#endif