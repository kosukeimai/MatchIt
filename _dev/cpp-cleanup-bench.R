#Benchmarks for the src/ cleanup. Run before and after; results go to
#_dev/cpp-cleanup-bench-<label>.rds. See _dev/cpp-cleanup-notes.md.
#PKG_DIR must have been compiled with the same optimization level on both sides:
#pkgbuild::compile_dll() defaults to -O0, which swamps everything measured here.
label <- Sys.getenv("BENCH_LABEL", "before")
pkg_dir <- Sys.getenv("PKG_DIR", "/Users/NoahGreifer/Dropbox/Research/R/MatchIt")

pkgload::load_all(pkg_dir, quiet = TRUE)
setwd("/Users/NoahGreifer/Dropbox/Research/R/MatchIt")

data("lalonde", package = "MatchIt")

set.seed(9876)
N <- 20000L
big <- data.frame(
  x1 = rnorm(N), x2 = rnorm(N), x3 = rbinom(N, 1, 0.4),
  x4 = sample(letters[1:4], N, TRUE), x5 = rpois(N, 3)
)
big$treat <- rbinom(N, 1, plogis(-1 + 0.5 * big$x1 - 0.3 * big$x2 + 0.4 * big$x3))
big$id <- rep(seq_len(N / 2L), each = 2L)

f_l <- treat ~ age + educ + race + married + nodegree + re74 + re75
f_b <- treat ~ x1 + x2 + x3 + x4 + x5

tm <- function(expr) {
  expr <- substitute(expr)
  min(replicate(3L, system.time(eval(expr, parent.frame()))[["elapsed"]]))
}

res <- c(
  #nn_matchC_vec: caliper + antiexact + exact
  nn_vec_lalonde = tm(matchit(f_l, data = lalonde, caliper = 0.2,
                              antiexact = ~race, exact = ~married)),
  nn_vec_big = tm(matchit(f_b, data = big, caliper = 0.2, antiexact = ~x4)),
  nn_vec_ratio = tm(matchit(f_b, data = big, ratio = 3, caliper = 0.5)),

  #nn_matchC_vec_closest / reuse.max
  nn_closest_big = tm(matchit(f_b, data = big, m.order = "closest",
                              caliper = 0.2)),
  nn_reuse_big = tm(matchit(f_b, data = big, replace = TRUE, reuse.max = 4,
                            ratio = 2)),

  #nn_matchC_mahcovs
  nn_mahcovs_lalonde = tm(matchit(f_l, data = lalonde,
                                  distance = "mahalanobis",
                                  antiexact = ~race)),
  nn_mahcovs_big = tm(matchit(f_b, data = big, distance = "mahalanobis")),
  nn_mahcovs_closest = tm(matchit(f_b, data = big, distance = "mahalanobis",
                                  m.order = "closest")),

  #unit.id path
  nn_unitid_big = tm(matchit(f_b, data = big, replace = TRUE, unit.id = ~id)),

  #optimal: preprocess_matchC
  optimal_lalonde = tm(matchit(f_l, data = lalonde, method = "optimal")),

  #stratification: subclass2mmC, weights_matrixC
  cem_big = tm(matchit(f_b, data = big, method = "cem")),
  full_lalonde = tm(matchit(f_l, data = lalonde, method = "full")),
  quick_big = tm(matchit(f_b, data = big, method = "quick")),

  #get_splitsC: caliper on a covariate with many unique values
  splits_big = tm(matchit(f_b, data = big, caliper = c(x1 = 0.05),
                          std.caliper = TRUE)),

  #pairdistsubC
  pairdist_full = tm(summary(matchit(f_l, data = lalonde, method = "full"),
                             pair.dist = TRUE))
)

print(round(res, 3))
saveRDS(res, sprintf("_dev/cpp-cleanup-bench-%s.rds", label))
