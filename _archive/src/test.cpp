// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

// [[Rcpp::plugins(cpp11)]]

struct nn_match2 : public Worker {
  const RVector<double> dist1;
  const RVector<double> dist0;
  const RVector<int> matched;
  const RVector<int> ind;
  const RVector<int> ind_un;

  RMatrix<int> pm;

  nn_match2(const Rcpp::NumericVector dist1,
            const Rcpp::NumericVector dist0,
            const Rcpp::IntegerVector matched,
            const Rcpp::IntegerVector ind,
            const Rcpp::IntegerVector ind_un,
            Rcpp::IntegerMatrix pm)
    : dist1(dist1), dist0(dist0), matched(matched), ind(ind),
      ind_un(ind_un), pm(pm) {}

  void operator()(std::size_t begin, std::size_t end) {
    int nu = ind.size();
    std::vector<double> distances(nu);
    int z, j, l, ind_l;

    for (int i = begin; i < end; i++) {
      std::vector<int> indices_sorted(i + 1);

      for (j = 0; j < nu; j++) {
        z = ind_un[j];
        distances[j] = std::abs(dist1[i] - dist0[z]);
      }

      std::partial_sort_copy(ind.begin(), ind.end(),
                             indices_sorted.begin(), indices_sorted.end(),
                             [&distances](int it0, int it1) {return distances[it0] < distances[it1];});

      for (l = 0; l < i + 1; l++) {
        ind_l = indices_sorted[l];
        pm(i, l) = ind_un[ind_l];
      }
    }
  }
};

// [[Rcpp::export]]
IntegerMatrix rcpp_parallel_nn_match2(const NumericVector& dist1,
                                      const NumericVector& dist0,
                                      const IntegerVector& matched) {

  // allocate the matrix we will return
  int n1 = dist1.size();
  int n0 = dist0.size();
  Rcpp::IntegerMatrix pm(n1, n0);
  IntegerVector ind_un = Range(0, n0 - 1);
  ind_un = ind_un[matched == 0];
  IntegerVector ind = Range(0, ind_un.size() - 1);

  // create the worker
  nn_match2 nnm2(dist1, dist0, matched, ind, ind_un, pm);

  // call it with parallelFor
  parallelFor(0, n1, nnm2, 10);

  return pm;
}

// [[Rcpp::export]]
IntegerMatrix nn_match2_w_par(const IntegerMatrix& mm_,
                              const NumericVector& dist1,
                              const NumericVector& dist0) {
  IntegerMatrix mm = mm_;

  int n1 = dist1.size();
  int n0 = dist0.size();
  IntegerVector matched = rep(0, n0);

  IntegerMatrix pm(n1, n0);

  pm = rcpp_parallel_nn_match2(dist1, dist0, matched);

  for (int i = 0; i < n1; i++) {
    IntegerVector i_ = pm.row(i);
    for (int j = 0; j < i + 1; j++) {
      int j_ = i_[j];
      if (matched[j_] == 0) {
        mm(i, 0) = j_;
        matched[j_] = 1;

        break;
      }
    }
  }

  return mm + 1;
}

// [[Rcpp::export]]
IntegerMatrix nn_match3_w_par(const IntegerMatrix& mm_,
                              const NumericVector& dist1,
                              const NumericVector& dist0,
                              const int& blocksize) {
  IntegerMatrix mm = mm_;

  int n1 = dist1.size();
  int n0 = dist0.size();
  IntegerVector matched = rep(0, n0);
  IntegerVector r(blocksize);
  IntegerMatrix pm;

  int done = 0;
  int k = 0;
  int newk;

  while (done == 0) {
    if (n1 - k > blocksize) newk = k + blocksize;
    else {
      newk = n1;
      done = 1;
    }

    r = Range(k, newk - 1);

    pm = rcpp_parallel_nn_match2(dist1[r], dist0, matched);

    for (int i = 0; i < r.size(); i++) {
      int ki = r[i];
      IntegerVector i_ = pm.row(i);
      for (int j = 0; j < i + 1; j++) {
        int j_ = i_[j];
        if (matched[j_] == 0) {
          mm(ki, 0) = j_;
          matched[j_] = 1;

          break;
        }
      }
    }

    k = newk;
  }

  return mm + 1;
}

// [[Rcpp::export]]
IntegerMatrix nn_match3(const IntegerMatrix& mm_,
                        const NumericVector& dist1,
                        const NumericVector& dist0,
                        const int& blocksize) {
  IntegerMatrix mm = mm_;

  int n1 = dist1.size();
  int n0 = dist0.size();
  IntegerVector matched = rep(0, n0);

  IntegerMatrix pm;
  NumericVector distances;

  int done = 0;
  int k = 0;
  int newk;
  int nu, i, j, ri, l, ind_l;

  while (done == 0) {
    if (n1 - k > blocksize) newk = k + blocksize;
    else {
      newk = n1;
      done = 1;
    }

    IntegerMatrix pm(newk - k, newk - k);

    nu = n0 - sum(matched);
    IntegerVector ind(nu), ind_un(nu);
    NumericVector distances(nu);
    IntegerVector indices_sorted;

    int z = 0;
    for (j = 0; j < n0; j++) {
      if (matched[j] == 0) {
        ind_un[z] = j;
        ind[z] = z;
        z++;
      }
    }

    for (i = 0; i < newk - k; i++) {
      ri = i + k;

      IntegerVector indices_sorted(i + 1);

      z = 0;
      for (j = 0; j < n0; j++) {
        if (matched[j] == 0) {
          distances[z] = abs(dist1[ri] - dist0[j]);
          z++;
        }
      }

      std::partial_sort_copy(ind.begin(), ind.end(),
                             indices_sorted.begin(), indices_sorted.end(),
                        [&distances](int it0, int it1) {return distances[it0] < distances[it1];});

      for (l = 0; l < i + 1; l++) {
        ind_l = indices_sorted[l];
        pm(i, l) = ind_un[ind_l];
      }
    }

    for (int i = 0; i < newk - k; i++) {
      ri = i + k;
      IntegerVector i_ = pm.row(i);
      for (int j = 0; j < i + 1; j++) {
        int j_ = i_[j];
        if (matched[j_] == 0) {
          mm(ri, 0) = j_;
          matched[j_] = 1;

          break;
        }
      }
    }

    k = newk;
  }

  return mm + 1;
}

// [[Rcpp::export]]
IntegerMatrix nn_match1(const IntegerMatrix& mm_,
                        const NumericVector& dist1,
                        const NumericVector& dist0) {
  IntegerMatrix mm = mm_;

  int n1 = dist1.size();
  int n0 = dist0.size();
  LogicalVector matched = rep(false, n0);

  IntegerVector indices(n0), eligible(n0);
  NumericVector distances(n0), edist(n0);
  int m, wm;

  //Order the matches up to the maximum needed to prevent repeats; parallelizable
  for (int i = 0; i < n1; i++) {
    indices = Range(0, n0 - 1);
    eligible = indices[!matched];
    edist = dist0[eligible];
    distances = abs(dist1[i] - edist);
    wm = which_min(distances);
    m = eligible[wm];
    mm(i,0) = m;
    matched[m] = true;
  }

  return mm + 1;
}

// [[Rcpp::export]]
IntegerMatrix nn_match2(const IntegerMatrix& mm_,
                        const NumericVector& dist1,
                        const NumericVector& dist0) {
  IntegerMatrix mm = mm_;

  int n1 = dist1.size();
  int n0 = dist0.size();
  LogicalVector matched = rep(false, n0);

  IntegerMatrix pm(n1, n0);
  IntegerVector indices(n0), r(n0);
  NumericVector distances(n0);

  //Order the matches up to the maximum needed to prevent repeats; parallelizable
  for (int i = 0; i < n1; i++) {
    distances = abs(dist1[i] - dist0);

    indices = Range(0, n0 - 1);

    std::partial_sort(indices.begin(), indices.begin() + i + 1, indices.end(),
                      [&distances](int k, int j) {return distances[k] < distances[j];});

    for (int l = 0; l < i + 1; l++) {
      pm(i, l) = indices[l];
    }
  }

  for (int i = 0; i < n1; i++) {
    IntegerVector i_ = pm.row(i);
    for (int j = 0; j < i + 1; j++) {
      int j_ = i_[j];
      if (!matched[j_]) {
        mm(i, 0) = j_;
        matched[j_] = true;

        break;
      }
    }
  }

  return mm + 1;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
d1 = rnorm(10000)
d0 = rnorm(20000)
m <- matrix(0L, nrow = length(d1), ncol = 1)

microbenchmark::microbenchmark(
  nn_match3(m,d1, d0, 100),
  # nn_match3(m,d1, d0, 50),
  # nn_match3(m,d1, d0, 25),
  # nn_match3(m,d1, d0, 10),
  # nn_match3(m,d1, d0, 5),
  nn_match3_w_par(m,d1, d0, 100),
  # nn_match3_w_par(m,d1, d0, 50),
  times = 10, check = "equivalent"
)

# do.call("rbind", lapply(seq_along(d1), \(i) {x = order(abs(d1[i] - d0))-1; x[-(1:i)] <- 0; x}))
*/
