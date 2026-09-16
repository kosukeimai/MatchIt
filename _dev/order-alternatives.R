#Are there faster replacements for the base::order() calls in src/, and do they give
#byte-identical results? Verifies first, then times. Companion to
#_dev/order-alternatives.cpp; conclusions in _dev/cpp-cleanup-notes.md.
#
#Run with: littler _dev/order-alternatives.R

Rcpp::sourceCpp("_dev/order-alternatives.cpp")

fns <- list(rcall = ord_rcall, capi = ord_capi, stdsort = ord_stdsort,
            stablesort = ord_stablesort, pairsort = ord_pairsort)

# ===== byte-identity against base::order() =====

set.seed(42)

cases <- list(
  continuous = rnorm(1000),
  few_unique = sample(c(0.1, 0.2, 0.3), 1000, TRUE),
  all_equal = rep(1.5, 500),
  two_values = rep(c(0, 1), 500),
  integers = as.numeric(sample(1:10, 1000, TRUE)),
  sorted = sort(rnorm(1000)),
  reverse = sort(rnorm(1000), decreasing = TRUE),
  ties_and_na = sample(c(1, 2, 2, 3, NA), 500, TRUE),
  na_nan = c(NA, NaN, 1, NaN, NA, 0.5, Inf, -Inf, 2),
  all_na = rep(NA_real_, 20),
  inf = c(Inf, -Inf, 0, Inf, -Inf, NA, NaN),
  neg_zero = c(0, -0, 1, -1, 0),
  tiny = c(1e-320, 0, -1e-320, 1e300),
  single = 1.5,
  empty = numeric(0),
  ps_like = round(plogis(rnorm(2000)), 6),
  discrete_ps = round(plogis(rnorm(2000)), 2)
)

bad <- 0L

check <- function(x) {
  for (dec in c(FALSE, TRUE)) {
    ref <- order(x, decreasing = dec) - 1L
    for (f in names(fns)) {
      if (!identical(as.integer(fns[[f]](x, dec)), as.integer(ref))) {
        bad <<- bad + 1L
        cat("MISMATCH:", f, "decreasing =", dec, "\n")
      }
    }
  }
}

for (x in cases) check(x)

#Fuzz, weighted toward ties and missingness
for (i in 1:300) {
  n <- sample(1:200, 1)
  check(switch(sample(3, 1),
               rnorm(n),
               sample(c(1, 2, 3, NA, NaN), n, TRUE),
               round(rnorm(n), 1)))
}

cat(if (bad == 0L) "All alternatives byte-identical to base::order()\n"
    else sprintf("%d mismatches\n", bad), "\n")

# ===== timing =====

bench1 <- function(x, reps) {
  vapply(fns, function(f) {
    min(replicate(5, system.time(for (i in seq_len(reps)) f(x))[["elapsed"]])) / reps * 1e6
  }, numeric(1))
}

sizes <- list(
  "n=614 continuous"     = list(x = rnorm(614), reps = 2000),
  "n=614 discrete(20)"   = list(x = as.numeric(sample(20, 614, TRUE)), reps = 2000),
  "n=20000 continuous"   = list(x = rnorm(20000), reps = 100),
  "n=20000 discrete(20)" = list(x = as.numeric(sample(20, 20000, TRUE)), reps = 100),
  "n=20000 presorted"    = list(x = sort(rnorm(20000)), reps = 100),
  "n=1e5 continuous"     = list(x = rnorm(1e5), reps = 30),
  "n=8000 heap-like"     = list(x = abs(rnorm(8000)), reps = 200)
)

cat("microseconds per call:\n")
print(round(t(vapply(sizes, function(s) bench1(s$x, s$reps), numeric(length(fns)))), 1))

#Why: order()'s default for a double vector is a radix sort; R_orderVector1() is the
#shell sort, which is what the C API exposes.
x <- rnorm(20000)
cat(sprintf("\norder(radix) %.0f us | order(shell) %.0f us\n",
            min(replicate(5, system.time(for (i in 1:100) order(x, method = "radix"))[["elapsed"]])) / 100 * 1e6,
            min(replicate(5, system.time(for (i in 1:100) order(x, method = "shell"))[["elapsed"]])) / 100 * 1e6))
