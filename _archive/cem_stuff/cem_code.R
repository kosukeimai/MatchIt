#R code that depends on cem. Return to indicated location when cem is returned to CRAN.

#Return to matchit2methods.R in section "MATCHIT method = cem"
matchit2cem <- function(treat, covs, estimand = "ATT", verbose = FALSE, ...) {

  if (length(covs) == 0) stop("Covariates must be specified in the input formula to use coarsened exact matching.", call. = FALSE)

  check.package("cem")

  if (verbose) cat("Coarsened exact matching...\n")

  A <- list(...)
  A[["method"]] <- A[["k2k.method"]]

  n.obs <- length(treat)

  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATT", "ATC", "ATE"))

  # cem takes the data all together and wants the treatment specified
  # with the column name of the data frame. Here we massage the matchit
  # inputs to this format. Note that X has its proper column names, but
  # treat does not have the original column name.
  cem.data <- data.frame(covs, treat)

  args.excluded <- c("treatment", "baseline.group", "data", "verbose", "eval.imbalance",
                     "keep.all", "drop", "L1.breaks", "L1.grouping")
  mat <- tryCatch({
    withCallingHandlers({
      do.call(cem::cem, c(list(treatment = names(cem.data)[ncol(cem.data)],
                               data = cem.data,
                               verbose = as.integer(verbose),
                               eval.imbalance = FALSE,
                               keep.all = FALSE,
                               drop = NULL),
                          A[names(A) %in% setdiff(names(formals(cem::cem)), args.excluded)]))
    },
    warning = function(w) {
      warning(paste0("(from cem) ", conditionMessage(w)), call. = FALSE, immediate. = TRUE)
      invokeRestart("muffleWarning")
    })
  },
  error = function(e) {
    if (startsWith(conditionMessage(e), "subscript out of bounds")) {
      stop("No units were matched. Try changing the coarsening options using the 'cutpoints' and 'grouping' arguments in cem(). See ?method_cem or ?cem::cem for details.", call. = FALSE)
    }
    else {
      stop(paste0("(from cem) ", conditionMessage(e)), call. = FALSE)
    }
  })

  strat <- setNames(rep(NA_character_, n.obs), names(treat))
  if (!is.null(mat)) strat[mat$matched] <- mat$strata[mat$matched]
  strat <- setNames(factor(strat, labels = seq_along(unique(strat[!is.na(strat)]))), names(treat))

  if (verbose) cat("Calculating matching weights... ")

  res <- list(subclass = strat,
              weights = weights.subclass(strat, treat, estimand))

  if (verbose) cat("Done.\n")

  class(res) <- "matchit"
  return(res)
}