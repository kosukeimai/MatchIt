#include <Rcpp.h>
#include <R_ext/Utils.h>
#include <algorithm>
#include <numeric>
using namespace Rcpp;

//1. What the package does now: call base::order() and shift to 0-based.
// [[Rcpp::export]]
IntegerVector ord_rcall(const NumericVector& x, bool decreasing = false) {
  Function o = Environment::base_env()["order"];
  IntegerVector out = o(x, _["decreasing"] = decreasing);
  return out - 1;
}

//2. R's own C entry point for order(), filling 0-based indices directly.
// [[Rcpp::export]]
IntegerVector ord_capi(const NumericVector& x, bool decreasing = false) {
  R_xlen_t n = x.size();
  IntegerVector out(n);
  R_orderVector1(out.begin(), n, x, TRUE, decreasing ? TRUE : FALSE);
  return out;
}

//3. std::sort over indices, with the index itself as the tiebreak so the result is
//   the stable ordering base::order() produces, and NaN/NA sorted last.
// [[Rcpp::export]]
IntegerVector ord_stdsort(const NumericVector& x, bool decreasing = false) {
  R_xlen_t n = x.size();
  IntegerVector out(n);
  std::iota(out.begin(), out.end(), 0);

  const double* px = x.begin();

  if (decreasing) {
    std::sort(out.begin(), out.end(), [px](int a, int b) {
      double xa = px[a], xb = px[b];
      bool na_a = ISNAN(xa), na_b = ISNAN(xb);
      if (na_a || na_b) {
        if (na_a && na_b) return a < b;
        return na_b;
      }
      if (xa > xb) return true;
      if (xb > xa) return false;
      return a < b;
    });
  }
  else {
    std::sort(out.begin(), out.end(), [px](int a, int b) {
      double xa = px[a], xb = px[b];
      bool na_a = ISNAN(xa), na_b = ISNAN(xb);
      if (na_a || na_b) {
        if (na_a && na_b) return a < b;
        return na_b;
      }
      if (xa < xb) return true;
      if (xb < xa) return false;
      return a < b;
    });
  }

  return out;
}

//4. stable_sort without the index tiebreak, for comparison.
// [[Rcpp::export]]
IntegerVector ord_stablesort(const NumericVector& x, bool decreasing = false) {
  R_xlen_t n = x.size();
  IntegerVector out(n);
  std::iota(out.begin(), out.end(), 0);

  const double* px = x.begin();

  std::stable_sort(out.begin(), out.end(), [px, decreasing](int a, int b) {
    double xa = px[a], xb = px[b];
    bool na_a = ISNAN(xa), na_b = ISNAN(xb);
    if (na_a || na_b) {
      if (na_a && na_b) return false;
      return na_b;
    }
    return decreasing ? (xa > xb) : (xa < xb);
  });

  return out;
}

//5. std::sort over (value, index) pairs: 16 bytes moved per element instead of an
//   indirection through x on every comparison.
// [[Rcpp::export]]
IntegerVector ord_pairsort(const NumericVector& x, bool decreasing = false) {
  R_xlen_t n = x.size();

  std::vector<std::pair<double, int>> v(n);
  for (R_xlen_t i = 0; i < n; i++) {
    v[i] = std::make_pair(x[i], static_cast<int>(i));
  }

  auto cmp = [decreasing](const std::pair<double, int>& a,
                          const std::pair<double, int>& b) {
    bool na_a = ISNAN(a.first), na_b = ISNAN(b.first);
    if (na_a || na_b) {
      if (na_a && na_b) return a.second < b.second;
      return na_b;
    }
    if (a.first != b.first) {
      return decreasing ? (a.first > b.first) : (a.first < b.first);
    }
    return a.second < b.second;
  };

  std::sort(v.begin(), v.end(), cmp);

  IntegerVector out(n);
  for (R_xlen_t i = 0; i < n; i++) {
    out[i] = v[i].second;
  }

  return out;
}
