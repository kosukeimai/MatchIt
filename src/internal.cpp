#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// Rcpp internal functions

//C implementation of tabulate(). Faster than base::tabulate(), but real
//use is in subclass2mmC().

// [[Rcpp::interfaces(cpp)]]
IntegerVector tabulateC_(const IntegerVector& bins,
                         const int& nbins = 0) {
  int max_bin;

  if (nbins > 0) max_bin = nbins;
  else max_bin = max(na_omit(bins));

  IntegerVector counts(max_bin);
  int n = bins.size();
  for (int i = 0; i < n; i++) {
    if (bins[i] > 0 && bins[i] <= max_bin) {
      counts[bins[i] - 1]++;
    }
  }

  return counts;
}

//Rcpp port of base::which

// [[Rcpp::interfaces(cpp)]]
IntegerVector which(const LogicalVector& x) {
  IntegerVector ind = Range(0, x.size() - 1);
  return ind[x];
}

// [[Rcpp::interfaces(cpp)]]
bool antiexact_okay(const int& aenc,
                    const int& i,
                    const int& j,
                    const IntegerMatrix& antiexact_covs) {
  if (aenc == 0) {
    return true;
  }

  IntegerVector antiexact_covs_row_i = antiexact_covs.row(i);
  IntegerVector antiexact_covs_row_j = antiexact_covs.row(j);

  for (int k = 0; k < aenc; k++) {
    if (antiexact_covs_row_i[k] == antiexact_covs(j, k)) {
      return false;
    }
  }

  return true;
}

// [[Rcpp::interfaces(cpp)]]
bool caliper_covs_okay(const int& ncc,
                       const int& i,
                       const int& j,
                       const NumericMatrix& caliper_covs_mat,
                       const NumericVector& caliper_covs) {
  if (ncc == 0) {
    return true;
  }

  for (int k = 0; k < ncc; k++) {
    if (caliper_covs[k] >= 0) {
      if (std::abs(caliper_covs_mat(i, k) - caliper_covs_mat(j, k)) > caliper_covs[k]) {
        return false;
      }
    }
    else {
      if (std::abs(caliper_covs_mat(i, k) - caliper_covs_mat(j, k)) <= -caliper_covs[k]) {
        return false;
      }
    }
  }

  return true;
}

// [[Rcpp::interfaces(cpp)]]
bool caliper_covs_okay2(const int& ncc,
                        const NumericVector& cc_ti,
                        const int& j,
                        const NumericMatrix& caliper_covs_mat,
                        const NumericVector& caliper_covs) {
  if (ncc == 0) {
    return true;
  }

  for (int k = 0; k < ncc; k++) {
    if (caliper_covs[k] >= 0) {
      if (std::abs(cc_ti[k] - caliper_covs_mat(j, k)) > caliper_covs[k]) {
        return false;
      }
    }
    else {
      if (std::abs(cc_ti[k] - caliper_covs_mat(j, k)) <= -caliper_covs[k]) {
        return false;
      }
    }
  }

  return true;
}

// [[Rcpp::interfaces(cpp)]]
bool caliper_dist_okay(const bool& use_caliper_dist,
                       const int& i,
                       const int& j,
                       const NumericVector& distance,
                       const double& caliper_dist) {
  if (!use_caliper_dist) {
    return true;
  }

  if (caliper_dist >= 0) {
    return std::abs(distance[i] - distance[j]) <= caliper_dist;
  }
  else {
    return std::abs(distance[i] - distance[j]) > -caliper_dist;
  }
}

// [[Rcpp::interfaces(cpp)]]
bool mm_okay(const int& r,
             const int& i,
             const IntegerVector& mm_rowi) {

  if (r > 1) {
    for (int j : na_omit(mm_rowi)) {
      if (i == j) {
        return false;
      }
    }
  }

  return true;
}

// [[Rcpp::interfaces(cpp)]]
bool exact_okay(const bool& use_exact,
                const int& i,
                const int& j,
                const IntegerVector& exact) {

  if (!use_exact) {
    return true;
  }

  return exact[i] == exact[j];
}

// [[Rcpp::interfaces(cpp)]]
double euc_dist_sq(const NumericVector& v1,
                   const NumericVector& v2) {
  double out = 0;
  double tmp = 0;
  int s = v1.size();

  for (int i = 0; i < s; i++) {
    tmp = v1[i] - v2[i];
    out += tmp * tmp;
  }

  return out;
}

// [[Rcpp::interfaces(cpp)]]
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
                                  const int& prev_start = -1) {

  int ii = match_d_ord[t_id];

  IntegerVector mm_rowi;
  std::vector<int> possible_starts;

  if (r > 1) {
    mm_rowi = na_omit(mm_rowi_);
    mm_rowi = mm_rowi[as<IntegerVector>(treat[mm_rowi]) == gi];
    possible_starts.reserve(mm_rowi.size() + 2);

    for (int mmi : mm_rowi) {
      possible_starts.push_back(match_d_ord[mmi]);
    }
  }
  else {
    possible_starts.reserve(2);
  }

  if (prev_start >= 0) {
    possible_starts.push_back(match_d_ord[prev_start]);
  }

  int iil, iir;
  double min_dist;

  if (possible_starts.empty()) {
    iil = ii;
    iir = ii;
    min_dist = 0;
  }
  else {
    possible_starts.push_back(ii);

    iil = *std::min_element(possible_starts.begin(), possible_starts.end());
    iir = *std::max_element(possible_starts.begin(), possible_starts.end());

    if (iil == ii) {
      min_dist = std::abs(distance[t_id] - distance[ind_d_ord[iir]]);
    }
    else if (iir == ii) {
      min_dist = std::abs(distance[t_id] - distance[ind_d_ord[iil]]);
    }
    else {
      min_dist = std::max(std::abs(distance[t_id] - distance[ind_d_ord[iil]]),
                          std::abs(distance[t_id] - distance[ind_d_ord[iir]]));
    }
  }

  if (caliper_dist <= 0 && min_dist < -caliper_dist) {
    min_dist = -caliper_dist;
  }

  int min_ii = first_control[gi];
  int max_ii = last_control[gi];

  double di = distance[t_id];

  NumericVector cc_ti;
  if (ncc > 0) {
    cc_ti = caliper_covs_mat.row(t_id);
  }

  bool l_stop = false;
  bool r_stop = false;

  double dist_c;

  std::vector<int> potential_matches_id;
  potential_matches_id.reserve(2 * ratio);
  std::vector<double> potential_matches_dist;
  potential_matches_dist.reserve(2 * ratio);

  int num_matches_l = 0;
  int num_matches_r = 0;

  int iz;
  bool left = false;
  int num_closer_than_dist_c;

  while (!l_stop || !r_stop) {
    if (l_stop) {
      left = false;
    }
    else if (r_stop) {
      left = true;
    }
    else {
      left = !left;
    }

    if (left) {
      if (iil <= min_ii || num_matches_l == ratio) {
        l_stop = true;
        continue;
      }

      iil -= 1;
      iz = ind_d_ord[iil];
    }
    else {
      if (iir >= max_ii || num_matches_r == ratio) {
        r_stop = true;
        continue;
      }

      iir += 1;
      iz = ind_d_ord[iir];
    }

    if (!eligible[iz]) {
      continue;
    }

    if (treat[iz] != gi) {
      continue;
    }

    if (!mm_okay(r, iz, mm_rowi)) {
      continue;
    }

    dist_c = std::abs(di - distance[iz]);

    if (caliper_dist >= 0) {
      if (dist_c > caliper_dist) {
        if (left) {
          l_stop = true;
        }
        else {
          r_stop = true;
        }
        continue;
      }
    }
    else if (dist_c <= -caliper_dist) {
      continue;
    }

    if (dist_c < min_dist) {
      continue;
    }

    //If current dist is worse than ratio dists, continue
    if (potential_matches_id.size() >= static_cast<size_t>(ratio)) {
      num_closer_than_dist_c = 0;
      for (double d : potential_matches_dist) {
        if (d < dist_c) {
          num_closer_than_dist_c++;
          if (num_closer_than_dist_c == ratio) {
            break;
          }
        }
      }

      if (num_closer_than_dist_c >= ratio) {
        if (left) {
          l_stop = true;
        }
        else {
          r_stop = true;
        }
        continue;
      }
    }



    if (!exact_okay(use_exact, t_id, iz, exact)) {
      continue;
    }

    if (!antiexact_okay(aenc, t_id, iz, antiexact_covs)) {
      continue;
    }

    // if (!caliper_covs_okay(ncc, t_id, iz, caliper_covs_mat, caliper_covs)) {
    //   continue;
    // }

    if (!caliper_covs_okay2(ncc, cc_ti, iz, caliper_covs_mat, caliper_covs)) {
      continue;
    }

    potential_matches_id.push_back(iz);
    potential_matches_dist.push_back(dist_c);

    if (left) {
      num_matches_l++;
      if (num_matches_l == ratio) {
        l_stop = true;
      }
    }
    else {
      num_matches_r++;
      if (num_matches_r == ratio) {
        r_stop = true;
      }
    }
  }

  int n_potential_matches = potential_matches_id.size();

  if (n_potential_matches <= 1) {
    return potential_matches_id;
  }

  if (n_potential_matches <= ratio &&
      std::is_sorted(potential_matches_dist.begin(),
                     potential_matches_dist.end())) {
    return potential_matches_id;
  }

  std::vector<int> ind(n_potential_matches);
  std::iota(ind.begin(), ind.end(), 0);

  std::vector<int> matches_out;

  if (n_potential_matches > ratio) {
    std::partial_sort(ind.begin(), ind.begin() + ratio, ind.end(),
                      [&potential_matches_dist](int a, int b){
                        return potential_matches_dist[a] < potential_matches_dist[b];
                      });

    matches_out.reserve(ratio);

    for (auto it = ind.begin(); it != ind.begin() + ratio; ++it) {
      matches_out.push_back(potential_matches_id[*it]);
    }
  }
  else {
    std::sort(ind.begin(), ind.end(),
              [&potential_matches_dist](int a, int b){
                return potential_matches_dist[a] < potential_matches_dist[b];
              });

    matches_out.reserve(n_potential_matches);

    for (auto it = ind.begin(); it != ind.end(); ++it) {
      matches_out.push_back(potential_matches_id[*it]);
    }
  }

  return matches_out;
}

// [[Rcpp::interfaces(cpp)]]
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
                                      const int& ratio = 1) {

  int ii = match_d_ord[t_id];

  int iil, iir;

  iil = ii;
  iir = ii;

  int min_ii = 0;
  int max_ii = match_d_ord.size() - 1;

  bool l_stop = false;
  bool r_stop = false;

  double dist_c;

  std::vector<std::pair<int, double>> potential_matches;
  potential_matches.reserve(ratio);

  std::pair<int,double> new_match;

  int num_matches_l = 0;
  int num_matches_r = 0;

  double mv_i = match_var[t_id];
  double mv_dist;

  int iz;
  bool left = false;

  auto dist_comp = [](std::pair<int, double> a, std::pair<int, double> b) {
    return a.second < b.second;
  };

  while (!l_stop || !r_stop) {
    if (l_stop) {
      left = false;
    }
    else if (r_stop) {
      left = true;
    }
    else {
      left = !left;
    }

    if (left) {
      if (iil <= min_ii || num_matches_l == ratio) {
        l_stop = true;
        continue;
      }

      iil -= 1;
      iz = ind_d_ord[iil];
    }
    else {
      if (iir >= max_ii || num_matches_r == ratio) {
        r_stop = true;
        continue;
      }

      iir += 1;
      iz = ind_d_ord[iir];
    }

    if (!eligible[iz]) {
      continue;
    }

    if (treat[iz] != gi) {
      continue;
    }

    if (!mm_okay(r, iz, mm_rowi)) {
      continue;
    }

    mv_dist = std::abs(mv_i - match_var[iz]);

    if (match_var_caliper >= 0) {
      if (mv_dist > match_var_caliper) {
        if (left) {
          l_stop = true;
        }
        else {
          r_stop = true;
        }
        continue;
      }
    }
    else if (mv_dist <= -match_var_caliper) {
      continue;
    }

    mv_dist = mv_dist * mv_dist;

    //If current dist is worse than ratio dists, continue
    if (potential_matches.size() == static_cast<size_t>(ratio)) {
      if (potential_matches.back().second < mv_dist) {
        if (left) {
          l_stop = true;
        }
        else {
          r_stop = true;
        }

        continue;
      }
    }

    if (!exact_okay(use_exact, t_id, iz, exact)) {
      continue;
    }

    if (!caliper_dist_okay(use_caliper_dist, t_id, iz, distance, caliper_dist)) {
      continue;
    }

    if (!antiexact_okay(aenc, t_id, iz, antiexact_covs)) {
      continue;
    }

    if (!caliper_covs_okay(ncc, t_id, iz, caliper_covs_mat, caliper_covs)) {
      continue;
    }

    dist_c = euc_dist_sq(mah_covs.row(t_id), mah_covs.row(iz));

    if (!std::isfinite(dist_c)) {
      continue;
    }

    new_match = std::pair<int,double>(iz, dist_c);

    if (potential_matches.empty()) {
      potential_matches.push_back(new_match);
    }
    else if (dist_c > potential_matches.back().second) {
      if (potential_matches.size() == static_cast<size_t>(ratio)) {
        continue;
      }

      potential_matches.push_back(new_match);
    }
    else if (ratio == 1) {
      potential_matches[0] = new_match;
    }
    else {
      if (potential_matches.size() == static_cast<size_t>(ratio)) {
        potential_matches.pop_back();
      }

      if (dist_c > potential_matches.back().second) {
        potential_matches.push_back(new_match);
      }
      else {
        potential_matches.insert(std::lower_bound(potential_matches.begin(), potential_matches.end(),
                                                  new_match, dist_comp),
                                                  new_match);
      }
    }
  }

  std::vector<int> matches_out;
  matches_out.reserve(potential_matches.size());

  for (auto p : potential_matches) {
    matches_out.push_back(p.first);
  }

  return matches_out;
}

// [[Rcpp::interfaces(cpp)]]
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
                                  const int& ratio = 1) {

  int c_id_i;
  double dist_c;

  std::vector<int> potential_matches_id;

  if (ratio < 1) {
    return potential_matches_id;
  }

  std::vector<double> potential_matches_dist;
  double max_dist = R_PosInf;

  R_xlen_t nc = distance_mat_row_i.size();

  potential_matches_id.reserve(nc);
  potential_matches_dist.reserve(nc);

  for (R_xlen_t c = 0; c < nc; c++) {

    dist_c = distance_mat_row_i[c];

    if (potential_matches_id.size() == static_cast<size_t>(ratio)) {
      if (dist_c > max_dist) {
        continue;
      }
    }

    if (caliper_dist >= 0) {
      if (dist_c > caliper_dist) {
        continue;
      }
    }
    else {
      if (dist_c <= -caliper_dist) {
        continue;
      }
    }

    if (!std::isfinite(dist_c)) {
      continue;
    }

    c_id_i = ind_non_focal[c];

    if (!eligible[c_id_i]) {
      continue;
    }

    if (treat[c_id_i] != gi) {
      continue;
    }

    if (!mm_okay(r, c_id_i, mm_rowi)) {
      continue;
    }

    if (!exact_okay(use_exact, t_id, c_id_i, exact)) {
      continue;
    }

    if (!antiexact_okay(aenc, t_id, c_id_i, antiexact_covs)) {
      continue;
    }

    if (!caliper_covs_okay(ncc, t_id, c_id_i, caliper_covs_mat, caliper_covs)) {
      continue;
    }

    potential_matches_id.push_back(c_id_i);
    potential_matches_dist.push_back(dist_c);

    if (potential_matches_id.size() == 1) {
      max_dist = dist_c;
    }
    else if (dist_c > max_dist) {
      max_dist = dist_c;
    }
  }

  int n_potential_matches = potential_matches_id.size();

  if (n_potential_matches <= 1) {
    return potential_matches_id;
  }

  if (n_potential_matches <= ratio &&
      std::is_sorted(potential_matches_dist.begin(),
                     potential_matches_dist.end())) {
    return potential_matches_id;
  }

  std::vector<int> ind(n_potential_matches);
  std::iota(ind.begin(), ind.end(), 0);

  std::vector<int> matches_out;

  if (n_potential_matches > ratio) {
    std::partial_sort(ind.begin(), ind.begin() + ratio, ind.end(),
                      [&potential_matches_dist](int a, int b){
                        return potential_matches_dist[a] < potential_matches_dist[b];
                      });

    matches_out.reserve(ratio);

    for (auto it = ind.begin(); it != ind.begin() + ratio; ++it) {
      matches_out.push_back(potential_matches_id[*it]);
    }
  }
  else {
    std::sort(ind.begin(), ind.end(),
              [&potential_matches_dist](int a, int b){
                return potential_matches_dist[a] < potential_matches_dist[b];
              });

    matches_out.reserve(n_potential_matches);

    for (auto it = ind.begin(); it != ind.end(); ++it) {
      matches_out.push_back(potential_matches_id[*it]);
    }
  }

  return matches_out;
}

// [[Rcpp::interfaces(cpp)]]
double max_finite(const NumericVector& x) {
  double m = NA_REAL;

  R_xlen_t n = x.size();
  R_xlen_t i;
  bool found = false;

  //Find first finite value
  for (i = 0; i < n; i++) {
    if (std::isfinite(x[i])) {
      m = x[i];
      found = true;
      break;
    }
  }

  //If none found, return NA
  if (!found) {
    return m;
  }

  //Find largest finite value
  for (R_xlen_t j = i + 1; j < n; j++) {
    if (!std::isfinite(x[j])) {
      continue;
    }

    if (x[j] > m) {
      m = x[j];
    }
  }

  return m;
}

// [[Rcpp::interfaces(cpp)]]
double min_finite(const NumericVector& x) {
  double m = NA_REAL;

  R_xlen_t n = x.size();
  R_xlen_t i;
  bool found = false;

  //Find first finite value
  for (i = 0; i < n; i++) {
    if (std::isfinite(x[i])) {
      m = x[i];
      found = true;
      break;
    }
  }

  //If none found, return NA
  if (!found) {
    return m;
  }

  //Find smallest finite value
  for (R_xlen_t j = i + 1; j < n; j++) {
    if (!std::isfinite(x[j])) {
      continue;
    }

    if (x[j] < m) {
      m = x[j];
    }
  }

  return m;
}

// [[Rcpp::interfaces(cpp)]]
void update_first_and_last_control(IntegerVector first_control,
                                   IntegerVector last_control,
                                   const IntegerVector& ind_d_ord,
                                   const LogicalVector& eligible,
                                   const IntegerVector& treat,
                                   const int& gi) {
  R_xlen_t c;

  // Update first_control
  if (!eligible[ind_d_ord[first_control[gi]]]) {
    for (c = first_control[gi] + 1; c <= last_control[gi]; c++) {
      if (eligible[ind_d_ord[c]]) {
        if (treat[ind_d_ord[c]] == gi) {
          first_control[gi] = c;
          break;
        }
      }
    }
  }

  // Update last_control
  if (!eligible[ind_d_ord[last_control[gi]]]) {
    for (c = last_control[gi] - 1; c >= first_control[gi]; c--) {
      if (eligible[ind_d_ord[c]]) {
        if (treat[ind_d_ord[c]] == gi) {
          last_control[gi] = c;
          break;
        }
      }
    }
  }
}

// [[Rcpp::interfaces(cpp)]]
double get_affine_transformation(const NumericVector& x,
                                 const NumericVector& y,
                                 const double& tol = 1e-9) {
  R_len_t n = x.size();
  int i;

  if (n != y.size() || n < 2) {
    return false; // Need at least two points for a meaningful check
  }

  // Compute means
  double mean_x = mean(x);
  double mean_y = mean(y);

  // Compute a (scaling factor)
  double num = 0.0, denom = 0.0;
  double x_diff, y_diff;
  for (i = 0; i < n; i++) {
    x_diff = x[i] - mean_x;
    y_diff = y[i] - mean_y;

    num += x_diff * y_diff;
    denom += x_diff * x_diff;
  }

  if (std::abs(denom) < tol || std::abs(num) < tol) {
    return 0.0;
  }

  double a = num / denom;
  double b = mean_y - a * mean_x;

  // Verify if y is reconstructed correctly within tolerance
  for (i = 0; i < n; i++) {
    if (std::abs(a * x[i] + b - y[i]) > tol) {
      return 0.0;
    }
  }

  return a;
}