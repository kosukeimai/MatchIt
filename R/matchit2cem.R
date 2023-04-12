#' Coarsened Exact Matching
#' @name method_cem
#' @aliases method_cem
#'
#' @usage NULL
#'
#' @description
#' In [matchit()], setting `method = "cem"` performs coarsened exact
#' matching. With coarsened exact matching, covariates are coarsened into bins,
#' and a complete cross of the coarsened covariates is used to form subclasses
#' defined by each combination of the coarsened covariate levels. Any subclass
#' that doesn't contain both treated and control units is discarded, leaving
#' only subclasses containing treatment and control units that are exactly
#' equal on the coarsened covariates. The coarsening process can be controlled
#' by an algorithm or by manually specifying cutpoints and groupings. The
#' benefits of coarsened exact matching are that the tradeoff between exact
#' matching and approximate balancing can be managed to prevent discarding too
#' many units, which can otherwise occur with exact matching.
#'
#' This page details the allowable arguments with `method = "cem"`. See
#' [matchit()] for an explanation of what each argument means in a general
#' context and how it can be specified.
#'
#' Below is how `matchit()` is used for coarsened exact matching:
#' \preformatted{
#' matchit(formula,
#'         data = NULL,
#'         method = "cem",
#'         estimand = "ATT",
#'         s.weights = NULL,
#'         verbose = FALSE,
#'         ...) }
#'
#' @param formula a two-sided [formula] object containing the treatment and
#' covariates to be used in creating the subclasses defined by a full cross of
#' the coarsened covariate levels.
#' @param data a data frame containing the variables named in `formula`.
#' If not found in `data`, the variables will be sought in the
#' environment.
#' @param method set here to `"cem"`.
#' @param estimand a string containing the desired estimand. Allowable options
#' include `"ATT"`, `"ATC"`, and `"ATE"`. The estimand controls
#' how the weights are computed; see the Computing Weights section at
#' [matchit()] for details. When `k2k = TRUE` (see below), `estimand`
#' also controls how the matching is done.
#' @param s.weights the variable containing sampling weights to be incorporated
#' into balance statistics or the scaling factors when `k2k = TRUE` and
#' certain methods are used.
#' @param verbose `logical`; whether information about the matching
#' process should be printed to the console.
#' @param \dots additional arguments to control the matching process.
#' \describe{
#' \item{`grouping`}{ a named list with an (optional) entry
#' for each categorical variable to be matched on. Each element should itself
#' be a list, and each entry of the sublist should be a vector containing
#' levels of the variable that should be combined to form a single level. Any
#' categorical variables not included in `grouping` will remain as they
#' are in the data, which means exact matching, with no coarsening, will take
#' place on these variables. See Details.  }
#' \item{`cutpoints`}{ a named list with an (optional) entry for each numeric variable to be matched on.
#' Each element describes a way of coarsening the corresponding variable. They
#' can be a vector of cutpoints that demarcate bins, a single number giving the
#' number of bins, or a string corresponding to a method of computing the
#' number of bins. Allowable strings include `"sturges"`, `"scott"`,
#' and `"fd"`, which use the functions
#' [grDevices::nclass.Sturges()], [grDevices::nclass.scott()],
#' and [grDevices::nclass.FD()], respectively. The default is
#' `"sturges"` for variables that are not listed or if no argument is
#' supplied. Can also be a single value to be applied to all numeric variables.
#' See Details.  }
#' \item{`k2k`}{ `logical`; whether 1:1 matching should
#' occur within the matched strata. If `TRUE` nearest neighbor matching
#' without replacement will take place within each stratum, and any unmatched
#' units will be dropped (e.g., if there are more treated than control units in
#' the stratum, the treated units without a match will be dropped). The
#' `k2k.method` argument controls how the distance between units is
#' calculated.  }
#' \item{`k2k.method`}{ `character`; how the distance
#' between units should be calculated if `k2k = TRUE`. Allowable arguments
#' include `NULL` (for random matching), any argument to
#' [distance()] for computing a distance matrix from covariates
#' (e.g., `"mahalanobis"`), or any allowable argument to `method` in
#' [dist()]. Matching will take place on the original
#' (non-coarsened) variables. The default is `"mahalanobis"`.  }
#' \item{`mpower`}{ if `k2k.method = "minkowski"`, the power used in
#' creating the distance. This is passed to the `p` argument of [dist()].}
#' }
#'
#' The arguments `distance` (and related arguments), `exact`, `mahvars`, `discard` (and related arguments), `replace`, `m.order`, `caliper` (and related arguments), and `ratio` are ignored with a warning.
#'
#' @section Outputs:
#'
#' All outputs described in [matchit()] are returned with
#' `method = "cem"` except for `match.matrix`. When `k2k = TRUE`, a `match.matrix` component with the matched pairs is also
#' included. `include.obj` is ignored.
#'
#' @details
#' If the coarsening is such that there are no exact matches with the coarsened
#' variables, the `grouping` and `cutpoints` arguments can be used to
#' modify the matching specification. Reducing the number of cutpoints or
#' grouping some variable values together can make it easier to find matches.
#' See Examples below. Removing variables can also help (but they will likely
#' not be balanced unless highly correlated with the included variables). To
#' take advantage of coarsened exact matching without failing to find any
#' matches, the covariates can be manually coarsened outside of
#' `matchit()` and then supplied to the `exact` argument in a call to
#' `matchit()` with another matching method.
#'
#' Setting `k2k = TRUE` is equivalent to first doing coarsened exact
#' matching with `k2k = FALSE` and then supplying stratum membership as an
#' exact matching variable (i.e., in `exact`) to another call to
#' `matchit()` with `method = "nearest"`, `distance =
#' "mahalanobis"` and an argument to `discard` denoting unmatched units.
#' It is also equivalent to performing nearest neighbor matching supplying
#' coarsened versions of the variables to `exact`, except that
#' `method = "cem"` automatically coarsens the continuous variables. The
#' `estimand` argument supplied with `method = "cem"` functions the
#' same way it would in these alternate matching calls, i.e., by determining
#' the "focal" group that controls the order of the matching.
#'
#' ## Grouping and Cutpoints
#'
#' The `grouping` and `cutpoints`
#' arguments allow one to fine-tune the coarsening of the covariates.
#' `grouping` is used for combining categories of categorical covariates
#' and `cutpoints` is used for binning numeric covariates. The values
#' supplied to these arguments should be iteratively changed until a matching
#' solution that balances covariate balance and remaining sample size is
#' obtained. The arguments are described below.
#'
#' ### `grouping`
#'
#' The argument to `grouping` must be a list, where each component has the
#' name of a categorical variable, the levels of which are to be combined. Each
#' component must itself be a list; this list contains one or more vectors of
#' levels, where each vector corresponds to the levels that should be combined
#' into a single category. For example, if a variable `amount` had levels
#' `"none"`, `"some"`, and `"a lot"`, one could enter
#' `grouping = list(amount = list(c("none"), c("some", "a lot")))`, which
#' would group `"some"` and `"a lot"` into a single category and
#' leave `"none"` in its own category. Any levels left out of the list for
#' each variable will be left alone (so `c("none")` could have been
#' omitted from the previous code). Note that if a categorical variable does
#' not appear in `grouping`, it will not be coarsened, so exact matching
#' will take place on it. `grouping` should not be used for numeric
#' variables with more than a few values; use `cutpoints`, described below, instead.
#'
#' ### `cutpoints`
#'
#' The argument to `cutpoints` must also be a list, where each component
#' has the name of a numeric variables that is to be binned. (As a shortcut, it
#' can also be a single value that will be applied to all numeric variables).
#' Each component can take one of three forms: a vector of cutpoints that
#' separate the bins, a single number giving the number of bins, or a string
#' corresponding to an algorithm used to compute the number of bins. Any values
#' at a boundary will be placed into the higher bin; e.g., if the cutpoints
#' were `c(0, 5, 10)`, values of 5 would be placed into the same bin as
#' values of 6, 7, 8, or 9, and values of 10 would be placed into a different
#' bin. Internally, values of `-Inf` and `Inf` are appended to the
#' beginning and end of the range. When given as a single number defining the
#' number of bins, the bin boundaries are the maximum and minimum values of the
#' variable with bin boundaries evenly spaced between them, i.e., not
#' quantiles. A value of 0 will not perform any binning (equivalent to exact
#' matching on the variable), and a value of 1 will remove the variable from
#' the exact matching variables but it will be still used for pair matching
#' when `k2k = TRUE`. The allowable strings include `"sturges"`,
#' `"scott"`, and `"fd"`, which use the corresponding binning method,
#' and `"q#"` where `#` is a number, which splits the variable into
#' `#` equally-sized bins (i.e., quantiles).
#'
#' An example of a way to supply an argument to `cutpoints` would be the
#' following:
#' \preformatted{
#' cutpoints = list(X1 = 4,
#'                  X2 = c(1.7, 5.5, 10.2),
#'                  X3 = "scott",
#'                  X4 = "q5") }
#'
#' This would split `X1` into 4 bins, `X2`
#' into bins based on the provided boundaries, `X3` into a number of bins
#' determined by [grDevices::nclass.scott()], and `X4` into
#' quintiles. All other numeric variables would be split into a number of bins
#' determined by [grDevices::nclass.Sturges()], the default.
#'
#'
#' @note
#' This method does not rely on the *cem* package, instead using
#' code written for *MatchIt*, but its design is based on the original
#' *cem* functions. Versions of *MatchIt* prior to 4.1.0 did rely on
#' *cem*, so results may differ between versions. There are a few
#' differences between the ways *MatchIt* and *cem* (and older
#' versions of *MatchIt*) differ in executing coarsened exact matching,
#' described below.
#' * In *MatchIt*, when a single number is
#' supplied to `cutpoints`, it describes the number of bins; in
#' *cem*, it describes the number of cutpoints separating bins. The
#' *MatchIt* method is closer to how [hist()] processes breaks points to
#' create bins.
#' * In *MatchIt*, values on the cutpoint boundaries will
#' be placed into the higher bin; in *cem*, they are placed into the lower
#' bin. To avoid consequences of this choice, ensure the bin boundaries do not
#' coincide with observed values of the variables.
#' * When `cutpoints` are used, `"ss"` (for Shimazaki-Shinomoto's rule) can be used in
#' *cem* but not in *MatchIt*.
#' * When `k2k = TRUE`, *MatchIt* matches on the original variables (scaled), whereas
#' *cem* matches on the coarsened variables. Because the variables are
#' already exactly matched on the coarsened variables, matching in *cem*
#' is equivalent to random matching within strata.
#' * When `k2k = TRUE`, in *MatchIt* matched units are identified by pair membership, and the
#' original stratum membership prior to 1:1 matching is discarded. In
#' *cem*, pairs are not identified beyond the stratum the members are part of.
#' * When `k2k = TRUE`, `k2k.method = "mahalanobis"` can be
#' requested in *MatchIt* but not in *cem*.
#'
#' @seealso [matchit()] for a detailed explanation of the inputs and outputs of
#' a call to `matchit()`.
#'
#' The *cem* package, upon which this method is based and which provided
#' the workhorse in previous versions of *MatchIt*.
#'
#' [`method_exact`] for exact matching, which performs exact matching
#' on the covariates without coarsening.
#'
#' @references
#' In a manuscript, you don't need to cite another package when
#' using `method = "cem"` because the matching is performed completely
#' within *MatchIt*. For example, a sentence might read:
#'
#' *Coarsened exact matching was performed using the MatchIt package (Ho,
#' Imai, King, & Stuart, 2011) in R.*
#'
#' It would be a good idea to cite the following article, which develops the
#' theory behind coarsened exact matching:
#'
#' Iacus, S. M., King, G., & Porro, G. (2012). Causal Inference without Balance
#' Checking: Coarsened Exact Matching. *Political Analysis*, 20(1), 1–24. \doi{10.1093/pan/mpr013}
#'
#' @examples
#'
#' data("lalonde")
#'
#' # Coarsened exact matching on age, race, married, and educ with educ
#' # coarsened into 5 bins and race coarsened into 2 categories,
#' # grouping "white" and "hispan" together
#' m.out1 <- matchit(treat ~ age + race + married + educ, data = lalonde,
#'                   method = "cem", cutpoints = list(educ = 5),
#'                   grouping = list(race = list(c("white", "hispan"),
#'                                               c("black"))))
#' m.out1
#' summary(m.out1)
#'
#' # The same but requesting 1:1 Mahalanobis distance matching with
#' # the k2k and k2k.method argument. Note the remaining number of units
#' # is smaller than when retaining the full matched sample.
#' m.out2 <- matchit(treat ~ age + race + married + educ, data = lalonde,
#'                   method = "cem", cutpoints = list(educ = 5),
#'                   grouping = list(race = list(c("white", "hispan"),
#'                                               "black")),
#'                   k2k = TRUE, k2k.method = "mahalanobis")
#' m.out2
#' summary(m.out2, un = FALSE)
#'
NULL

matchit2cem <- function(treat, covs, estimand = "ATT", s.weights = NULL, verbose = FALSE, ...) {
  if (length(covs) == 0) {
    .err("Covariates must be specified in the input formula to use coarsened exact matching")
  }

  if (verbose) cat("Coarsened exact matching...\n")

  A <- list(...)

  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATT", "ATC", "ATE"))

  #Uses in-house cem, no need for cem package. See cem_matchit.R for code.
  strat <- do.call("cem_matchit", c(list(treat = treat, X = covs, estimand = estimand,
                                         s.weights = s.weights),
                                    A[names(A) %in% names(formals(cem_matchit))]),
                   quote = TRUE)

  levels(strat) <- seq_len(nlevels(strat))
  names(strat) <- names(treat)

  mm <- NULL
  if (isTRUE(A[["k2k"]])) {
    mm <- nummm2charmm(subclass2mmC(strat, treat, focal = switch(estimand, "ATC" = 0, 1)),
                       treat)
  }

  if (verbose) cat("Calculating matching weights... ")

  res <- list(match.matrix = mm,
              subclass = strat,
              weights = weights.subclass(strat, treat, estimand))

  if (verbose) cat("Done.\n")

  class(res) <- "matchit"

  res
}

cem_matchit <- function(treat, X, cutpoints = "sturges", grouping = list(), k2k = FALSE,
                        k2k.method = "mahalanobis", mpower = 2, s.weights = NULL,
                        estimand = "ATT") {
  #In-house implementation of cem. Basically the same except:
  #treat is a vector if treatment status, not the name of a variable
  #X is a data.frame of covariates
  #when cutpoints are given as integer or string, they define the number of bins, not the number of breakpoints. "ss" is no longer allowed.
  #When k2k = TRUE, subclasses are created for each pair, mimicking true matching, not each covariate combo.
  #k2k.method is used instead of method. When k2k.method = NULL, units are matched based on order rather than random. Default is "mahalanobis" (not available in cem).
  #k2k now works with single covariates (previously it was ignored). k2k uses original variables, not coarsened versions

  if (k2k) {
    if (length(unique(treat)) > 2) {
      .err("`k2k` cannot be `TRUE` with a multi-category treatment")
    }
    if (!is.null(k2k.method)) {
    k2k.method <- tolower(k2k.method)
    k2k.method <- match_arg(k2k.method, c(matchit_distances(), "maximum", "manhattan", "canberra", "binary", "minkowski"))

    X.match <- transform_covariates(data = X, s.weights = s.weights, treat = treat,
                                    method = if (k2k.method %in% matchit_distances()) k2k.method else "euclidean")
    }
  }

  for (i in names(X)) {
    if (is.ordered(X[[i]])) X[[i]] <- as.numeric(X[[i]])
  }
  is.numeric.cov <- setNames(vapply(X, is.numeric, logical(1L)), names(X))

  #Process grouping
  if (length(grouping) > 0) {
    if (!is.list(grouping) || is.null(names(grouping))) {
      .err("`grouping` must be a named list of grouping values with an element for each variable whose values are to be binned")
    }

    bad.names <- setdiff(names(grouping), names(X))
    nb <- length(bad.names)
    if (nb > 0) {
      .wrn(sprintf("the variable%%s %s named in `grouping` %%r not in the variables supplied to `matchit()` and will be ignored",
                   word_list(bad.names, quotes = 2, and.or = "and")), n = nb)
      grouping[bad.names] <- NULL
    }

    for (i in names(grouping)) {
      X[[i]] <- as.character(X[[i]])
    }

    bag.groupings <- names(grouping)[vapply(grouping, function(g) {
      !is.list(g) ||
        !all(vapply(g, function(gg) is.atomic(gg) && is.vector(gg), logical(1L)))
    }, logical(1L))]
    nbg <- length(bag.groupings)
    if (nbg > 0) {
      .err(paste0("Each entry in the list supplied to `groupings` must be a list with entries containing values of the corresponding variable.",
                  "\nIncorrectly specified variable%s:\n\t"),
           paste(bag.groupings, collapse = ", "),
           tidy = FALSE, n = nbg)
    }

    for (g in names(grouping)) {
      x <- X[[g]]
      groups <- grouping[[g]]

      for (i in seq_along(groups)) {
        groups[[i]] <- as.character(groups[[i]])
        x[x %in% groups[[i]]] <- groups[[i]][1]
      }
      X[[g]] <- x

      #Remove cutpoints if variable named in `grouping`
      is.numeric.cov[g] <- FALSE
    }
  }

  #Process cutpoints
  if (!is.list(cutpoints)) {
    cutpoints <- setNames(lapply(names(X)[is.numeric.cov], function(i) cutpoints), names(X)[is.numeric.cov])
  }

  if (is.null(names(cutpoints))) {
    .err("`cutpoints` must be a named list of binning values with an element for each numeric variable")
  }
  bad.names <- setdiff(names(cutpoints), names(X))
  nb <- length(bad.names)
  if (nb > 0) {
    .wrn(sprintf("the variable%%s %s named in `cutpoints` %%r not in the variables supplied to `matchit()` and will be ignored",
                 word_list(bad.names, quotes = 2, and.or = "and")), n = nb)
    cutpoints[bad.names] <- NULL
  }

  if (length(grouping) > 0) {
    grouping.cutpoint.names <- intersect(names(grouping), names(cutpoints))
    ngc <- length(grouping.cutpoint.names)
    if (ngc > 0) {
      .wrn(sprintf("the variable%%s %s %%r named in both `grouping` and `cutpoints`; %s entr%%y%%s in `cutpoints` will be ignored",
                   word_list(grouping.cutpoint.names, quotes = 2, and.or = "and"),
                   ngettext(ngc, "its", "their")), n = ngc)
      cutpoints[grouping.cutpoint.names] <- NULL
    }
  }

  non.numeric.in.cutpoints <- intersect(names(X)[!is.numeric.cov], names(cutpoints))
  nnnic <- length(non.numeric.in.cutpoints)
  if (nnnic > 0) {
    .wrn(sprintf("the variable%%s %s named in `cutpoints` %%r not numeric and %s cutpoints will not be applied. Use `grouping` for non-numeric variables",
                 word_list(non.numeric.in.cutpoints, quotes = 2, and.or = "and"),
                 ngettext(nnnic, "its", "their")), n = nnnic)
  }

  bad.cuts <- setNames(rep(FALSE, length(cutpoints)), names(cutpoints))
  for (i in names(cutpoints)) {
    if (length(cutpoints[[i]]) == 0) {
      cutpoints[[i]] <- "sturges"
    }
    else if (length(cutpoints[[i]]) == 1) {
      if (is.na(cutpoints[[i]])) is.numeric.cov[i] <- FALSE #Will not be binned
      else if (is.character(cutpoints[[i]])) {
        bad.cuts[i] <- !(startsWith(cutpoints[[i]], "q") && can_str2num(substring(cutpoints[[i]], 2))) &&
          is.na(pmatch(cutpoints[[i]], c("sturges", "fd", "scott")))
      }
      else if (is.numeric(cutpoints[[i]])) {
        if (!is.finite(cutpoints[[i]]) || cutpoints[[i]] < 0) {
          bad.cuts[i] <- TRUE
        }
        else if (cutpoints[[i]] == 0) {
          is.numeric.cov[i] <- FALSE #Will not be binned
        }
        else if (cutpoints[[i]] == 1) {
          X[[i]] <- NULL #Removing from X, still in X.match
          is.numeric.cov <- is.numeric.cov[names(is.numeric.cov) != i]
        }
      }
      else {
        bad.cuts[i] <- TRUE
      }
    }
    else {
      bad.cuts[i] <- !is.numeric(cutpoints[[i]])
    }
  }

  if (any(bad.cuts)) {
    .err(paste0("All entries in the list supplied to `cutpoints` must be one of the following:",
                "\n\t- a string containing the name of an allowable binning method",
                "\n\t- a single number corresponding to the number of bins",
                "\n\t- a numeric vector containing the cut points separating bins",
                "\nIncorrectly specified variable%s:\n\t"),
                paste(names(cutpoints)[bad.cuts], collapse = ", "),
         tidy = FALSE, n = sum(bad.cuts))
  }

  #Create bins for numeric variables
  for (i in names(X)[is.numeric.cov]) {
    if (is.null(cutpoints) || !i %in% names(cutpoints)) bins <- "sturges"
    else bins <- cutpoints[[i]]

    if (is.character(bins)) {
      if (startsWith(bins, "q") || can_str2num(substring(bins, 2))) {
        #Quantile bins
        q <- str2num(substring(bins, 2))
        bins <- quantile(X[[i]], probs = seq(1/q, 1 - 1/q, by = 1/q), names = FALSE) #Outer boundaries will be added later
      }
      else {
        bins <- match_arg(tolower(bins), c("sturges", "fd", "scott"))
        bins <- switch(bins,
                       sturges = nclass.Sturges(X[[i]]),
                       fd = nclass.FD(X[[i]]),
                       scott = nclass.scott(X[[i]]))
        #Breaks is now a single number
      }
    }

    if (length(bins) == 1) {
      #cutpoints is number of bins, unlike in cem
      breaks <- seq(min(X[[i]]), max(X[[i]]), length = bins + 1)
      breaks[c(1, bins + 1)] <- c(-Inf, Inf)
    }
    else {
      breaks <- c(-Inf, sort(unique(bins)), Inf)
    }

    X[[i]] <- findInterval(X[[i]], breaks)
  }

  if (length(X) == 0) {
    subclass <- setNames(rep(1, length(treat)), names(treat))
  }
  else {
    #Exact match
    xx <- exactify(X, names(treat))
    cc <- do.call("intersect", unname(split(xx, treat)))

    if (length(cc) == 0) {
      .err("no units were matched. Try coarsening the variables further or decrease the number of variables to match on")
    }

    subclass <- setNames(match(xx, cc), names(treat))
  }

  extra.sub <- max(subclass, na.rm = TRUE)

  if (k2k) {

    na.sub <- is.na(subclass)

    s <- switch(estimand, "ATC" = 0, 1)

    for (i in which(tabulateC(subclass[!na.sub]) > 2)) {

      in.sub <- which(!na.sub & subclass == i)

      #Compute distance matrix; all 0s if k2k.method = NULL for matched based on data order
      if (is.null(k2k.method)) {
        dist.mat <- matrix(0, nrow = sum(treat[in.sub] == s), ncol = sum(treat[in.sub] != s),
                           dimnames = list(names(treat)[in.sub][treat[in.sub] == s],
                                           names(treat)[in.sub][treat[in.sub] != s]))

      }
      else if (k2k.method %in% matchit_distances()) {
        #X.match has been transformed
        dist.mat <- eucdist_internal(X.match[in.sub,,drop = FALSE], treat[in.sub] == s)
      }
      else {
        dist.mat <- dist_to_matrixC(dist(X.match[in.sub,,drop = FALSE], method = k2k.method, p = mpower))

        #Put smaller group on rows
        d.rows <- which(rownames(dist.mat) %in% names(treat[in.sub])[treat[in.sub] == s])
        dist.mat <- dist.mat[d.rows, -d.rows, drop = FALSE]
      }

      #For each member of group on row, find closest remaining pair from cols
      while (all(dim(dist.mat) > 0)) {
        extra.sub <- extra.sub + 1

        closest <- which.min(dist.mat[1,])
        subclass[c(rownames(dist.mat)[1], colnames(dist.mat)[closest])] <- extra.sub

        #Drop already paired units from dist.mat
        dist.mat <- dist.mat[-1,-closest, drop = FALSE]
      }

      #If any unmatched units remain, give them NA subclass
      if (any(dim(dist.mat) > 0)) is.na(subclass)[unlist(dimnames(dist.mat))] <- TRUE

    }
  }

  subclass <- factor(subclass, nmax = extra.sub)

  setNames(subclass, names(treat))
}