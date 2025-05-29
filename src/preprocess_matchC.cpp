#include "internal.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

//Preprocess by pruning unnecessary edges as in Sävje (2020) https://doi.org/10.1214/19-STS699.
//Returns a vector of matrix indices for n1xn0 distance matrix.

// [[Rcpp::export]]
IntegerVector preprocess_matchC(IntegerVector t,
                                NumericVector p) {
  R_xlen_t n = t.size();
  R_xlen_t n1 = std::count(t.begin(), t.end(), 1);
  R_xlen_t n0 = n - n1;

  R_xlen_t i, j;

  Function ord("order");

  IntegerVector o = ord(p);
  o = o - 1; //location of each unit after sorting

  IntegerVector im(n);
  int i0 = 0;
  int i1 = 0;
  for (j = 0; j < n; j++) {
    if (t[j] == 0) {
      im[j] = i0;
      i0++;
    }
    else {
      im[j] = i1;
      i1++;
    }
  }

  IntegerVector a(n0), b(n0);

  std::vector<int> queue;
  queue.reserve(n1);

  int k0 = 0;
  int ci = 0;

  for (i = 0; i < n; i++) {
    if (t[o[i]] == 1) {
      queue.push_back(i);
      continue;
    }

    if (queue.size() > static_cast<size_t>(k0)) {
      a[ci] = queue[k0];
      k0++;
    }
    else {
      a[ci] = i;
    }

    ci++;
  }

  k0 = 0;
  ci = n0 - 1;

  queue.clear();

  for (i = n - 1; i >= 0; i--) {
    if (t[o[i]] == 1) {
      queue.push_back(i);
      continue;
    }

    if (queue.size() > static_cast<size_t>(k0)) {
      b[ci] = queue[k0];
      k0++;
    }
    else {
      b[ci] = i;
    }

    ci--;
  }

  std::vector<int> keep;
  keep.reserve(n1 * n0);

  ci = 0;

  for (i = 0; i < n; i++) {
    if (t[o[i]] == 1) {
      continue;
    }

    for (j = a[ci]; j <= b[ci]; j++) {
      if (t[o[j]] == 0) {
        continue;
      }

      keep.push_back(im[o[j]] + n1 * im[o[i]]);
    }

    ci++;
  }

  int keep_size = keep.size();
  IntegerVector out(keep_size);

  for (i = 0; i < keep_size; i++) {
    out[i] = keep[i] + 1;
  }

  return out;
}
