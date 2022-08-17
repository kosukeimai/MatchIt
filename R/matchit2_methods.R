# this function takes inputs from matchit() and returns the
# strata for each observation in the subclass entry, the
# weight for each observation in the weight entry, and the
# match.matrix object
#
# MATCHIT method = NULL--------------------------------------
matchit2null <- function(discarded, ...) {

  res <- list(weights = as.numeric(!discarded))
  class(res) <- "matchit"

  return(res)
}
# MATCHIT method = cem--------------------------------------
matchit2cem <- function(treat, covs, estimand = "ATT", s.weights = NULL, verbose = FALSE, ...) {
  if (length(covs) == 0) stop("Covariates must be specified in the input formula to use coarsened exact matching.", call. = FALSE)

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
  return(res)
}
# MATCHIT method = exact------------------------------------
matchit2exact <- function(treat, covs, data, estimand = "ATT", verbose = FALSE, ...){

  if(verbose)
    cat("Exact matching... \n")

  if (length(covs) == 0) stop("Covariates must be specified in the input formula to use exact matching.", call. = FALSE)

  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATT", "ATC", "ATE"))

  xx <- exactify(covs, names(treat))
  cc <- do.call("intersect", lapply(unique(treat), function(t) xx[treat == t]))

  if (length(cc) == 0) {
    stop("No exact matches were found.", call. = FALSE)
  }

  psclass <- setNames(factor(match(xx, cc), nmax = length(cc)), names(treat))

  if (verbose) cat("Calculating matching weights... ")

  res <- list(subclass = psclass,
              weights = weights.subclass(psclass, treat, estimand))

  if (verbose) cat("Done.\n")

  class(res) <- "matchit"
  return(res)
}

# MATCHIT method = full-------------------------------------
matchit2full <- function(treat, formula, data, distance, discarded,
                         ratio = NULL, s.weights = NULL, #min.controls and max.controls in attrs of replace
                         caliper = NULL, mahvars = NULL, exact = NULL,
                         estimand = "ATT", verbose = FALSE,
                         is.full.mahalanobis, antiexact = NULL, ...) {

  check.package("optmatch")

  if (verbose) cat("Full matching... \n")

  A <- list(...)

  fm.args <- c("omit.fraction", "mean.controls", "tol", "solver")
  A[!names(A) %in% fm.args] <- NULL

  #Set max problem size to Inf and return to original value after match
  omps <- getOption("optmatch_max_problem_size")
  on.exit(options(optmatch_max_problem_size = omps))
  options(optmatch_max_problem_size = Inf)

  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATT", "ATC", "ATE"))
  if (estimand == "ATC") {
    tc <- c("control", "treated")
    focal <- 0
  }
  else {
    tc <- c("treated", "control")
    focal <- 1
  }

  treat_ <- setNames(as.integer(treat[!discarded] == focal), names(treat)[!discarded])

  # if (!is.null(data)) data <- data[!discarded,]

  if (is.full.mahalanobis) {
    if (length(attr(terms(formula, data = data), "term.labels")) == 0) {
      stop(sprintf("Covariates must be specified in the input formula when distance = \"%s\".",
                   attr(is.full.mahalanobis, "transform")), call. = FALSE)
    }
    mahvars <- formula
  }

  min.controls <- attr(ratio, "min.controls")
  max.controls <- attr(ratio, "max.controls")

  #Exact matching strata
  if (!is.null(exact)) {
    ex <- factor(exactify(model.frame(exact, data = data),
                          sep = ", ", include_vars = TRUE)[!discarded])
    cc <- intersect(ex[treat_==1], ex[treat_==0])
    if (length(cc) == 0) stop("No matches were found.", call. = FALSE)
  }
  else {
    ex <- factor(rep("_", length(treat_)), levels = "_")
  }

  #Create distance matrix; note that Mahalanobis distance computed using entire
  #sample (minus discarded), like method2nearest, as opposed to within exact strata, like optmatch.
  if (!is.null(mahvars)) {
    transform <- if (is.full.mahalanobis) attr(is.full.mahalanobis, "transform") else "mahalanobis"
    mahcovs <- transform_covariates(mahvars, data = data, method = transform,
                                    s.weights = s.weights, treat = treat,
                                    discarded = discarded)
    mo <- eucdist_internal(mahcovs, treat)
  }
  else if (is.matrix(distance)) {
    mo <- distance
  }
  else {
    mo <- eucdist_internal(setNames(distance, names(treat)), treat)
  }

  #Transpose distance mat as needed
  if (focal == 0) mo <- t(mo)

  #Remove discarded units from distance mat
  mo <- mo[!discarded[treat == focal], !discarded[treat != focal], drop = FALSE]
  dimnames(mo) <- list(names(treat_)[treat_ == 1], names(treat_)[treat_ == 0])

  mo <- optmatch::match_on(mo, data = data[!discarded,, drop = FALSE])
  mo <- optmatch::as.InfinitySparseMatrix(mo)

  #Process antiexact
  if (!is.null(antiexact)) {
    antiexactcovs <- model.frame(antiexact, data)
    for (i in seq_len(ncol(antiexactcovs))) {
      mo <- mo + optmatch::antiExactMatch(antiexactcovs[[i]][!discarded], z = treat_)
    }
  }

  #Process caliper
  if (!is.null(caliper)) {
    if (min.controls != 0) {
      stop("Calipers cannot be used with method = \"full\" when 'min.controls' is specified.", call. = FALSE)
    }

    if (any(names(caliper) != "")) {
      cov.cals <- setdiff(names(caliper), "")
      calcovs <- get.covs.matrix(reformulate(cov.cals, intercept = FALSE), data = data)
    }
    for (i in seq_along(caliper)) {
      if (names(caliper)[i] != "") {
        mo_cal <- optmatch::match_on(setNames(calcovs[!discarded, names(caliper)[i]], names(treat_)), z = treat_)
      }
      else if (is.null(mahvars) || is.matrix(distance)) {
        mo_cal <- mo
      }
      else {
        mo_cal <- optmatch::match_on(setNames(distance[!discarded], names(treat_)), z = treat_)
      }

      mo <- mo + optmatch::caliper(mo_cal, caliper[i])
    }
    rm(mo_cal)
  }

  #Initialize pair membership; must include names
  pair <- setNames(rep(NA_character_, length(treat)), names(treat))
  p <- setNames(vector("list", nlevels(ex)), levels(ex))

  t_df <- data.frame(treat)

  for (e in levels(ex)) {
    if (nlevels(ex) > 1) {
      mo_ <- mo[ex[treat_==1] == e, ex[treat_==0] == e]
    }
    else mo_ <- mo

    if (any(dim(mo_) == 0) || !any(is.finite(mo_))) next
    else if (all(dim(mo_) == 1)) {
      pair[ex == e] <- paste(1, e, sep = "|")
      next
    }

    withCallingHandlers({
      p[[e]] <- do.call(optmatch::fullmatch,
                        c(list(mo_,
                               min.controls = min.controls,
                               max.controls = max.controls,
                               data = t_df), #just to get rownames; not actually used in matching
                          A))
    },
    warning = function(w) {
      warning(paste0("(from optmatch) ", conditionMessage(w)), call. = FALSE, immediate. = TRUE)
      invokeRestart("muffleWarning")
    },
    error = function(e1) {
      stop(paste0("(from optmatch) ", conditionMessage(e1)), call. = FALSE)
    })

    pair[names(p[[e]])[!is.na(p[[e]])]] <- paste(as.character(p[[e]][!is.na(p[[e]])]), e, sep = "|")
  }

  if (all(is.na(pair))) stop("No matches were found.", call. = FALSE)
  if (length(p) == 1) p <- p[[1]]

  psclass <- factor(pair)
  levels(psclass) <- seq_len(nlevels(psclass))
  names(psclass) <- names(treat)

  #No match.matrix because treated units don't index matched strata (i.e., more than one
  #treated unit can be in the same stratum). Stratum information is contained in subclass.

  if (verbose) cat("Calculating matching weights... ")

  res <- list(subclass = psclass,
              weights = weights.subclass(psclass, treat, estimand),
              obj = p)

  if (verbose) cat("Done.\n")

  class(res) <- c("matchit")
  return(res)
}

# MATCHIT method = optimal----------------------------------
matchit2optimal <- function(treat, formula, data, distance, discarded,
                            ratio = 1, s.weights = NULL, caliper = NULL,
                            mahvars = NULL, exact = NULL,
                            estimand = "ATT", verbose = FALSE,
                            is.full.mahalanobis,  antiexact = NULL, ...) {

  check.package("optmatch")

  if (verbose) cat("Optimal matching... \n")

  A <- list(...)
  pm.args <- c("tol", "solver")
  A[!names(A) %in% pm.args] <- NULL

  #Set max problem size to Inf and return to original value after match
  omps <- getOption("optmatch_max_problem_size")
  on.exit(options(optmatch_max_problem_size = omps))
  options(optmatch_max_problem_size = Inf)

  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATT", "ATC"))
  if (estimand == "ATC") {
    tc <- c("control", "treated")
    focal <- 0
  }
  else {
    tc <- c("treated", "control")
    focal <- 1
  }

  treat_ <- setNames(as.integer(treat[!discarded] == focal), names(treat)[!discarded])

  # if (!is.null(data)) data <- data[!discarded,]

  if (is.full.mahalanobis) {
    if (length(attr(terms(formula, data = data), "term.labels")) == 0) {
      stop(sprintf("Covariates must be specified in the input formula when distance = \"%s\".",
                   attr(is.full.mahalanobis, "transform")), call. = FALSE)
    }
    mahvars <- formula
  }

  if (!is.null(caliper)) {
    warning("Calipers are currently not compatible with method = \"optimal\" and will be ignored.", call. = FALSE, immediate. = TRUE)
    caliper <- NULL
  }

  min.controls <- attr(ratio, "min.controls")
  max.controls <- attr(ratio, "max.controls")

  if (is.null(max.controls)) {
    min.controls <- max.controls <- ratio
  }

  #Exact matching strata
  if (!is.null(exact)) {
    ex <- factor(exactify(model.frame(exact, data = data),
                          sep = ", ", include_vars = TRUE)[!discarded])

    cc <- intersect(ex[treat_==1], ex[treat_==0])
    if (length(cc) == 0) stop("No matches were found.", call. = FALSE)

    e_ratios <- vapply(levels(ex), function(e) sum(treat_[ex == e] == 0)/sum(treat_[ex == e] == 1), numeric(1L))

    if (any(e_ratios < 1)) {
      warning(sprintf("Fewer %s units than %s units in some 'exact' strata; not all %s units will get a match.",
                      tc[2], tc[1], tc[1]), immediate. = TRUE, call. = FALSE)
    }
    if (ratio > 1 && any(e_ratios < ratio)) {
      if (ratio == max.controls)
        warning(sprintf("Not all %s units will get %s matches.",
                        tc[1], ratio), immediate. = TRUE, call. = FALSE)
      else
        warning(sprintf("Not enough %s units for an average of %s matches per %s unit in all 'exact' strata.",
                        tc[2], ratio, tc[1]), immediate. = TRUE, call. = FALSE)
    }
  }
  else {
    ex <- factor(rep("_", length(treat_)), levels = "_")
    e_ratios <- setNames(sum(treat_ == 0)/sum(treat_ == 1), levels(ex))

    if (e_ratios < 1) {
      warning(sprintf("Fewer %s units than %s units; not all %s units will get a match.",
                      tc[2], tc[1], tc[1]), immediate. = TRUE, call. = FALSE)
    }
    else if (e_ratios < ratio) {
      if (ratio == max.controls)
        warning(sprintf("Not all %s units will get %s matches.",
                        tc[1], ratio), immediate. = TRUE, call. = FALSE)
      else
        warning(sprintf("Not enough %s units for an average of %s matches per %s unit.",
                        tc[2], ratio, tc[1]), immediate. = TRUE, call. = FALSE)
    }
  }

  #Create distance matrix; note that Mahalanobis distance computed using entire
  #sample (minus discarded), like method2nearest, as opposed to within exact strata, like optmatch.
  if (!is.null(mahvars)) {
    transform <- if (is.full.mahalanobis) attr(is.full.mahalanobis, "transform") else "mahalanobis"
    mahcovs <- transform_covariates(mahvars, data = data, method = transform,
                                    s.weights = s.weights, treat = treat,
                                    discarded = discarded)
    mo <- eucdist_internal(mahcovs, treat)
  }
  else if (is.matrix(distance)) {
    mo <- distance
  }
  else {
    mo <- eucdist_internal(setNames(distance, names(treat)), treat)
  }

  #Transpose distance mat as needed
  if (focal == 0) mo <- t(mo)

  #Remove discarded units from distance mat
  mo <- mo[!discarded[treat == focal], !discarded[treat != focal], drop = FALSE]
  dimnames(mo) <- list(names(treat_)[treat_ == 1], names(treat_)[treat_ == 0])

  mo <- optmatch::match_on(mo, data = data[!discarded,, drop = FALSE])
  mo <- optmatch::as.InfinitySparseMatrix(mo)

  #Process antiexact
  if (!is.null(antiexact)) {
    antiexactcovs <- model.frame(antiexact, data)
    for (i in seq_len(ncol(antiexactcovs))) {
      mo <- mo + optmatch::antiExactMatch(antiexactcovs[[i]][!discarded], z = treat_)
    }
  }

  #Initialize pair membership; must include names
  pair <- setNames(rep(NA_character_, length(treat)), names(treat))
  p <- setNames(vector("list", nlevels(ex)), levels(ex))

  t_df <- data.frame(treat)

  for (e in levels(ex)) {
    if (nlevels(ex) > 1) {
      mo_ <- mo[ex[treat_==1] == e, ex[treat_==0] == e]
    }
    else mo_ <- mo

    if (any(dim(mo_) == 0)) next
    else if (all(dim(mo_) == 1)) {
      pair[ex == e] <- paste(1, e, sep = "|")
      next
    }

    #Process ratio, etc., when available ratio in exact matching categories
    #(e_ratio) differs from requested ratio
    if (e_ratios[e] < 1) {
      #Switch treatment and control labels; unmatched treated units are dropped
      ratio_ <- min.controls_ <- max.controls_ <- 1
      mo_ <- t(mo_)
    }
    else if (e_ratios[e] < ratio) {
      #Lower ratio and min.controls.
      ratio_ <- e_ratios[e]
      min.controls_ <- min(min.controls, floor(e_ratios[e]))
      max.controls_ <- max.controls
    }
    else {
      ratio_ <- ratio
      min.controls_ <- min.controls
      max.controls_ <- max.controls
    }

    withCallingHandlers({
      p[[e]] <- do.call(optmatch::fullmatch,
                        c(list(mo_,
                               mean.controls = ratio_,
                               min.controls = min.controls_,
                               max.controls = max.controls_,
                               data = t_df), #just to get rownames; not actually used in matching
                          A))
    },
    warning = function(w) {
      warning(paste0("(from optmatch) ", conditionMessage(w)), call. = FALSE, immediate. = TRUE)
      invokeRestart("muffleWarning")
    },
    error = function(e1) {
      stop(paste0("(from optmatch) ", conditionMessage(e1)), call. = FALSE)
    })

    pair[names(p[[e]])[!is.na(p[[e]])]] <- paste(as.character(p[[e]][!is.na(p[[e]])]), e, sep = "|")
  }

  if (all(is.na(pair))) stop("No matches were found.", call. = FALSE)
  if (length(p) == 1) p <- p[[1]]

  psclass <- factor(pair)
  levels(psclass) <- seq_len(nlevels(psclass))
  names(psclass) <- names(treat)

  mm <- nummm2charmm(subclass2mmC(psclass, treat, focal), treat)

  if (verbose) cat("Calculating matching weights... ")

  ## calculate weights and return the results
  res <- list(match.matrix = mm,
              subclass = psclass,
              weights = weights.subclass(psclass, treat, estimand),
              obj = p)

  if (verbose) cat("Done.\n")

  class(res) <- "matchit"
  return(res)
}
# MATCHIT method = genetic----------------------------------
matchit2genetic <- function(treat, data, distance, discarded,
                            ratio = 1, s.weights = NULL, replace = FALSE, m.order = NULL,
                            caliper = NULL, mahvars = NULL, exact = NULL,
                            formula = NULL, estimand = "ATT", verbose = FALSE,
                            is.full.mahalanobis, use.genetic = TRUE,
                            antiexact = NULL, ...) {

  check.package(c("Matching", "rgenoud"))

  if (verbose) cat("Genetic matching... \n")

  A <- list(...)

  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATT", "ATC"))
  if (estimand == "ATC") {
    tc <- c("control", "treated")
    focal <- 0
  }
  else {
    tc <- c("treated", "control")
    focal <- 1
  }

  if (!replace) {
    if (sum(!discarded & treat != focal) < sum(!discarded & treat == focal)) {
      warning(sprintf("Fewer %s units than %s units; not all %s units will get a match.",
                      tc[2], tc[1], tc[1]), immediate. = TRUE, call. = FALSE)
    }
    else if (sum(!discarded & treat != focal) < sum(!discarded & treat == focal)*ratio) {
      stop(sprintf("Not enough %s units for %s matches for each %s unit.",
                   tc[2], ratio, tc[1]), call. = FALSE)
    }
  }

  treat <- setNames(as.integer(treat == focal), names(treat))

  n.obs <- length(treat)
  n1 <- sum(treat == 1)

  if (is.null(names(treat))) names(treat) <- seq_len(n.obs)

  m.order <- {
    if (is.null(distance)) match_arg(m.order, c("data", "random"))
    else if (!is.null(m.order)) match_arg(m.order, c("largest", "smallest", "random", "data"))
    else if (estimand == "ATC") "smallest"
    else "largest"
  }

  ord <- switch(m.order,
                "largest" = order(distance, decreasing = TRUE),
                "smallest" = order(distance),
                "random" = sample.int(n.obs),
                "data" = seq_len(n.obs))
  ord <- ord[!ord %in% which(discarded)]

  #Create X (matching variables) and covs_to_balance
  covs_to_balance <- get.covs.matrix(formula, data = data)
  if (!is.null(mahvars)) {
    X <- get.covs.matrix.for.dist(mahvars, data = data)
  }
  else if (is.full.mahalanobis) {
    X <- covs_to_balance
  }
  else {
    X <- cbind(covs_to_balance, distance)
  }

  if (ncol(covs_to_balance) == 0) {
    stop("Covariates must be specified in the input formula to use genetic matching.", call. = FALSE)
  }

  #Process exact; exact.log will be supplied to GenMatch() and Match()
  if (!is.null(exact)) {
    #Add covariates in exact not in X to X
    ex <- as.integer(factor(exactify(model.frame(exact, data = data), names(treat), sep = ", ", include_vars = TRUE)))

    cc <- intersect(ex[treat==1], ex[treat==0])
    if (length(cc) == 0) stop("No matches were found.", call. = FALSE)

    X <- cbind(X, ex)

    exact.log <- c(rep(FALSE, ncol(X) - 1), TRUE)
  }
  else exact.log <- ex <- NULL

  #Process caliper; cal will be supplied to GenMatch() and Match()
  if (!is.null(caliper)) {
    #Add covariates in caliper other than distance (cov.cals) not in X to X
    cov.cals <- setdiff(names(caliper), "")
    if (length(cov.cals) > 0 && any(!cov.cals %in% colnames(X))) {
      calcovs <- get.covs.matrix(reformulate(cov.cals[!cov.cals %in% colnames(X)]), data = data)
      X <- cbind(X, calcovs)

      #Expand exact.log for newly added covariates
      if (!is.null(exact.log)) exact.log <- c(exact.log, rep(FALSE, ncol(calcovs)))
    }
    else cov.cals <- NULL

    #Matching::Match multiplies calipers by pop SD, so we need to divide by pop SD to unstandardize
    pop.sd <- function(x) sqrt(sum((x-mean(x))^2)/length(x))
    caliper <- caliper / vapply(names(caliper), function(x) {
      if (x == "") pop.sd(distance[!discarded])
      else pop.sd(X[!discarded, x])
    }, numeric(1L))

    #cal needs one value per variable in X
    cal <- setNames(rep(Inf, ncol(X)), colnames(X))

    #First put covariate calipers into cal
    if (length(cov.cals) > 0) {
      cal[intersect(cov.cals, names(cal))] <- caliper[intersect(cov.cals, names(cal))]
    }

    #Then put distance caliper into cal
    if ("" %in% names(caliper)) {
      dist.cal <- caliper[names(caliper) == ""]
      if (!is.null(mahvars)) {
        #If mahvars specified, distance is not yet in X, so add it to X
        X <- cbind(X, distance)
        cal <- c(cal, dist.cal)

        #Expand exact.log for newly added distance
        if (!is.null(exact.log)) exact.log <- c(exact.log, FALSE)
      }
      else {
        #Otherwise, distance is in X at the specified index
        cal[ncol(covs_to_balance) + 1] <- dist.cal
      }
    }
    else dist.cal <- NULL

  }
  else {
    cal <- dist.cal <- cov.cals <- NULL
  }

  #Reorder data according to m.order since Match matches in order of data;
  #ord already excludes discarded units

  treat_ <- treat[ord]
  covs_to_balance <- covs_to_balance[ord,,drop = FALSE]
  X <- X[ord,,drop = FALSE]
  if (!is.null(s.weights)) s.weights <- s.weights[ord]

  if (!is.null(antiexact)) {
    antiexactcovs <- model.frame(antiexact, data)[ord,,drop = FALSE]
    antiexact_restrict <- cbind(do.call("rbind", lapply(seq_len(ncol(antiexactcovs)), function(i) {
      unique.vals <- unique(antiexactcovs[,i])
      do.call("rbind", lapply(unique.vals, function(u) {
        t(combn(which(antiexactcovs[,i] == u), 2))
      }))
    })), -1)
    if (!is.null(A[["restrict"]])) A[["restrict"]] <- rbind(A[["restrict"]], antiexact_restrict)
    else A[["restrict"]] <- antiexact_restrict
  }
  else {
    antiexactcovs <- NULL
  }

  if (use.genetic) {
    withCallingHandlers({
      g.out <- do.call(Matching::GenMatch,
                       c(list(Tr = treat_, X = X, BalanceMatrix = covs_to_balance,
                              M = ratio, exact = exact.log, caliper = cal,
                              replace = replace, estimand = "ATT", ties = FALSE,
                              CommonSupport = FALSE, verbose = verbose,
                              weights = s.weights, print.level = 2*verbose),
                         A[names(A) %in% names(formals(Matching::GenMatch))]))
    },
    warning = function(w) {
      if (!startsWith(conditionMessage(w), "replace==FALSE, but there are more (weighted) treated obs than control obs.")) {
        warning(paste0("(from Matching) ", conditionMessage(w)), call. = FALSE, immediate. = TRUE)
      }
      invokeRestart("muffleWarning")
    },
    error = function(e) {
      stop(paste0("(from Matching) ", conditionMessage(e)), call. = FALSE)
    })
  }
  else {
    #For debugging
    g.out <- NULL
  }

  lab <- names(treat)
  lab1 <- names(treat[treat == 1])

  lab_ <- names(treat_)

  ind_ <- seq_along(treat)[ord]

  # if (!isFALSE(A$use.Match)) {
  withCallingHandlers({
    m.out <- Matching::Match(Tr = treat_, X = X,
                             M = ratio, exact = exact.log, caliper = cal,
                             replace = replace, estimand = "ATT", ties = FALSE,
                             weights = s.weights, CommonSupport = FALSE, Weight = 3,
                             Weight.matrix = if (use.genetic) g.out
                             else if (is.null(s.weights)) generalized_inverse(cor(X))
                             else generalized_inverse(cov.wt(X, s.weights, cor = TRUE)$cor),
                             restrict = A[["restrict"]], version = "fast")
  },
  warning = function(w) {
    if (!startsWith(conditionMessage(w), "replace==FALSE, but there are more (weighted) treated obs than control obs.")) {
      warning(paste0("(from Matching) ", conditionMessage(w)), call. = FALSE, immediate. = TRUE)
    }
    invokeRestart("muffleWarning")
  },
  error = function(e) {
    stop(paste0("(from Matching) ", conditionMessage(e)), call. = FALSE)
  })

  #Note: must use character match.matrix because of re-ordering treat into treat_
  mm <- matrix(NA_integer_, nrow = n1, ncol = max(table(m.out$index.treated)),
               dimnames = list(lab1, NULL))

  unique.matched.focal <- unique(m.out$index.treated, nmax = n1)

  ind1__ <- match(lab_, lab1)
  for (i in unique.matched.focal) {
    matched.units <- ind_[m.out$index.control[m.out$index.treated == i]]
    mm[ind1__[i], seq_along(matched.units)] <- matched.units
  }

  # }
  # else {
  #   #Use nn_match() instead of Match()
  #   ord1 <- ord[ord %in% which(treat == 1)]
  #   if (!is.null(cov.cals)) calcovs <- get.covs.matrix(reformulate(cov.cals), data = data)
  #   else calcovs <- NULL
  #
  #   if (is.null(g.out)) MWM <- generalized_inverse(cov(X))
  #   else MWM <- g.out$Weight.matrix %*% diag(1/apply(X, 2, var))
  #
  #   if (isFALSE(A$fast)) {
  #     mm <- nn_match(treat, ord1, ratio, replace, discarded, distance, ex, dist.cal,
  #                    cov.cals, calcovs, X, MWM)
  #   }
  #   else {
  #     mm <- nn_matchC(treat, ord1, ratio, replace, discarded, distance, ex, dist.cal,
  #                     cov.cals, calcovs, X, MWM)
  #   }
  #
  #   mm[] <- names(treat)[mm]
  #   dimnames(mm) <- list(lab1, seq_len(ratio))
  # }

  if (verbose) cat("Calculating matching weights... ")

  if (replace) {
    psclass <- NULL
  }
  else {
    psclass <- mm2subclass(mm, treat)
  }

  res <- list(match.matrix = nummm2charmm(mm, treat),
              subclass = psclass,
              weights = weights.matrix(mm, treat),
              obj = g.out)

  if (verbose) cat("Done.\n")

  class(res) <- "matchit"
  return(res)
}

# MATCHIT method = nearest----------------------------------
matchit2nearest <-  function(treat, data, distance, discarded,
                             ratio = 1, s.weights = NULL, replace = FALSE, m.order = NULL,
                             caliper = NULL, mahvars = NULL, exact = NULL,
                             formula = NULL, estimand = "ATT", verbose = FALSE,
                             is.full.mahalanobis, fast = TRUE,
                             antiexact = NULL, unit.id = NULL, ...){

  if (verbose) {
    if (fast) check.package("RcppProgress")
    cat("Nearest neighbor matching... \n")
  }

  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATT", "ATC"))
  if (estimand == "ATC") {
    tc <- c("control", "treated")
    focal <- 0
  }
  else {
    tc <- c("treated", "control")
    focal <- 1
  }

  treat <- setNames(as.integer(treat == focal), names(treat))

  if (is.full.mahalanobis) {
    if (length(attr(terms(formula, data = data), "term.labels")) == 0) {
      stop(sprintf("Covariates must be specified in the input formula when distance = \"%s\".",
                   attr(is.full.mahalanobis, "transform")), call. = FALSE)
    }
    mahvars <- formula
  }

  n.obs <- length(treat)
  n1 <- sum(treat == 1)
  n0 <- n.obs - n1

  lab <- names(treat)
  lab1 <- lab[treat == 1]

  if (!is.null(distance)) {
    names(distance) <- names(treat)
  }

  min.controls <- attr(ratio, "min.controls")
  max.controls <- attr(ratio, "max.controls")

  mahcovs <- mahSigma_inv <- distance_mat <- NULL
  if (!is.null(mahvars)) {
    transform <- if (is.full.mahalanobis) attr(is.full.mahalanobis, "transform") else "mahalanobis"
    mahcovs <- transform_covariates(mahvars, data = data, method = transform,
                                    s.weights = s.weights, treat = treat,
                                    discarded = discarded)
  }
  else if (is.matrix(distance)) {
    distance_mat <- distance
    distance <- NULL
  }

  #Process caliper
  if (!is.null(caliper)) {
    if (any(names(caliper) != "")) {
      caliper.covs <- caliper[names(caliper) != ""]
      caliper.covs.mat <- get.covs.matrix(reformulate(names(caliper.covs)), data = data)
    }
    else {
      caliper.covs.mat <- caliper.covs <- NULL
    }

    if (any(names(caliper) == "")) {
      caliper.dist <- caliper[names(caliper) == ""]
    }
    else {
      caliper.dist <- NULL
    }
  }
  else {
    caliper.dist <- caliper.covs <- NULL
    caliper.covs.mat <- NULL
  }

  #Process antiexact
  if (!is.null(antiexact)) {
    antiexactcovs <- model.frame(antiexact, data)
    antiexactcovs <- do.call("cbind", lapply(seq_len(ncol(antiexactcovs)), function(i) {
      as.integer(as.factor(antiexactcovs[[i]]))
    }))
  }
  else {
    antiexactcovs <- NULL
  }

  reuse.max <- attr(replace, "reuse.max")

  if (reuse.max >= n1) {
    m.order <- "data"
  }

  if (!is.null(unit.id) && reuse.max < n1) {
    unit.id <- process.variable.input(unit.id, data)
    unit.id <- factor(exactify(model.frame(unit.id, data = data),
                             nam = lab, sep = ", ", include_vars = TRUE))
    num_ctrl_unit.ids <- length(unique(unit.id[treat == 0]))

    #If each control unit is a unit.id, unit.ids are meaningless
    if (num_ctrl_unit.ids == n0) unit.id <- NULL
  }
  else {
    unit.id <- NULL
  }

  if (!is.null(exact)) {
    ex <- factor(exactify(model.frame(exact, data = data), nam = lab, sep = ", ", include_vars = TRUE))

    cc <- intersect(as.integer(ex)[treat==1], as.integer(ex)[treat==0])
    if (length(cc) == 0) stop("No matches were found.", call. = FALSE)

    if (reuse.max < n1) {

      e_ratios <- vapply(levels(ex), function(e) {
        if (is.null(unit.id)) sum(treat[ex == e] == 0)*(reuse.max/sum(treat[ex == e] == 1))
        else length(unique(unit.id[treat == 0 & ex == e]))*(reuse.max/sum(treat[ex == e] == 1))
      }, numeric(1L))

      if (any(e_ratios < 1)) {
        warning(sprintf("Fewer %s units than %s units in some 'exact' strata; not all %s units will get a match.",
                       tc[2], tc[1], tc[1]), immediate. = TRUE, call. = FALSE)
      }
      if (ratio > 1 && any(e_ratios < ratio)) {
        if (is.null(max.controls) || ratio == max.controls)
          warning(sprintf("Not all %s units will get %s matches.",
                          tc[1], ratio), immediate. = TRUE, call. = FALSE)
        else
          warning(sprintf("Not enough %s units for an average of %s matches per %s unit in all 'exact' strata.",
                          tc[2], ratio, tc[1]), immediate. = TRUE, call. = FALSE)
      }
    }
  }
  else {
    ex <- NULL

    if (reuse.max < n1) {

      e_ratios <- {
        if (is.null(unit.id)) as.numeric(reuse.max)*n0/n1
        else as.numeric(reuse.max)*num_ctrl_unit.ids/n1
      }

      if (e_ratios < 1) {
        warning(sprintf("Fewer %s %s than %s units; not all %s units will get a match.",
                        tc[2], if (is.null(unit.id)) "units" else "unit IDs", tc[1], tc[1]),
                immediate. = TRUE, call. = FALSE)
      }
      else if (e_ratios < ratio) {
        if (is.null(max.controls) || ratio == max.controls)
          warning(sprintf("Not all %s units will get %s matches.",
                          tc[1], ratio),
                  immediate. = TRUE, call. = FALSE)
        else
          warning(sprintf("Not enough %s %s for an average of %s matches per %s unit.",
                          tc[2], if (is.null(unit.id)) "units" else "unit IDs", ratio, tc[1]),
                  immediate. = TRUE, call. = FALSE)
      }
    }
  }

  #Variable ratio (extremal matching), Ming & Rosenbaum (2000)
  #Each treated unit get its own value of ratio
  if (!is.null(max.controls)) {
    if (is.null(distance)) {
      if (is.full.mahalanobis) stop(sprintf("'distance' cannot be \"%s\" for variable ratio matching.",
                                            transform), call. = FALSE)
      else stop("'distance' cannot be supplied as a matrix for variable ratio matching.", call. = FALSE)
    }

    m <- round(ratio * n1)
    # if (m > sum(treat == 0)) stop("'ratio' must be less than or equal to n0/n1.", call. = FALSE)

    kmax <- floor((m - min.controls*(n1-1)) / (max.controls - min.controls))
    kmin <- n1 - kmax - 1
    kmed <- m - (min.controls*kmin + max.controls*kmax)

    ratio0 <- c(rep(min.controls, kmin), kmed, rep(max.controls, kmax))

    #Make sure no units are assigned 0 matches
    if (any(ratio0 == 0)) {
      ind <- which(ratio0 == 0)
      ratio0[ind] <- 1
      ratio0[ind + 1] <- ratio0[ind + 1] - 1
    }

    ratio <- rep(NA_integer_, n1)

    #Order by distance; treated are supposed to have higher values
    ratio[order(distance[treat == 1],
                decreasing = mean(distance[treat == 1]) > mean(distance[treat != 1]))] <- ratio0
    ratio <- as.integer(ratio)
  }
  else {
    ratio <- as.integer(rep(ratio, n1))
  }
  max_rat <- max(ratio)

  m.order <- {
    if (is.null(distance)) match_arg(m.order, c("data", "random"))
    else if (!is.null(m.order)) match_arg(m.order, c("largest", "smallest", "random", "data"))
    else if (estimand == "ATC") "smallest"
    else "largest"
  }

  if (is.null(ex) || !is.null(unit.id)) {
    ord <- switch(m.order,
                  "largest" = order(distance[treat == 1], decreasing = TRUE),
                  "smallest" = order(distance[treat == 1], decreasing = FALSE),
                  "random" = sample.int(n1),
                  "data" = seq_len(n1))
    mm <- nn_matchC(treat, ord, ratio, max_rat, discarded, reuse.max, distance, distance_mat, ex, caliper.dist,
                    caliper.covs, caliper.covs.mat, mahcovs, antiexactcovs, unit.id, verbose)
  }
  else {
    distance_ <- caliper.covs.mat_ <- mahcovs_ <- antiexactcovs_ <- distance_mat_ <- NULL
    mm_list <- lapply(levels(ex), function(e) {
      if (verbose) {
        cat(sprintf("Matching subgroup %s/%s: %s...\n",
                    match(e, levels(ex)), nlevels(ex), e))
      }

      .e <- which(ex == e)
      treat_ <- treat[.e]
      discarded_ <- discarded[.e]
      if (!is.null(distance)) distance_ <- distance[.e]
      if (!is.null(caliper.covs.mat)) caliper.covs.mat_ <- caliper.covs.mat[.e,,drop = FALSE]
      if (!is.null(mahcovs)) mahcovs_ <- mahcovs[.e,,drop = FALSE]
      if (!is.null(antiexactcovs)) antiexactcovs_ <- antiexactcovs[.e,,drop = FALSE]
      if (!is.null(distance_mat)) {
        .e1 <- which(ex[treat==1] == e)
        .e0 <- which(ex[treat==0] == e)
        distance_mat_ <- distance_mat[.e1, .e0, drop = FALSE]
      }
      ratio_ <- ratio[ex[treat==1]==e]

      n1_ <- sum(treat_ == 1)
      ord_ <- switch(m.order,
                     "largest" = order(distance_[treat_ == 1], decreasing = TRUE),
                     "smallest" = order(distance_[treat_ == 1], decreasing = FALSE),
                     "random" = sample.int(n1_),
                     "data" = seq_len(n1_))

      mm_ <- nn_matchC(treat_, ord_, ratio_, max_rat, discarded_, reuse.max, distance_, distance_mat_, NULL, caliper.dist,
                       caliper.covs, caliper.covs.mat_, mahcovs_, antiexactcovs_, NULL, verbose)

      #Ensure matched indices correspond to indices in full sample, not subgroup
      mm_[] <- seq_along(treat)[.e][mm_]
      mm_
    })

    mm <- do.call("rbind", mm_list)[lab1,, drop = FALSE]
  }

  if (verbose) cat("Calculating matching weights... ")

  if (reuse.max > 1) {
    psclass <- NULL
  }
  else {
    psclass <- mm2subclass(mm, treat)
  }

  res <- list(match.matrix = nummm2charmm(mm, treat),
              subclass = psclass,
              weights = weights.matrix(mm, treat))

  if (verbose) cat("Done.\n")

  class(res) <- "matchit"

  return(res)
}

# MATCHIT method = subclass---------------------------------
matchit2subclass <- function(treat, distance, discarded,
                             replace = FALSE, exact = NULL,
                             estimand = "ATT", verbose = FALSE,
                             ...) {

  if(verbose)
    cat("Subclassifying... \n")

  A <- list(...)
  subclass <- A[["subclass"]]
  sub.by <- A[["sub.by"]]
  min.n <- A[["min.n"]]

  #Checks
  if (is.null(subclass)) subclass <- 6
  else if (!is.numeric(subclass) || !is.null(dim(subclass))) {
    stop("subclass must be a numeric value.", call. = FALSE)
  }
  else if (length(subclass) == 1) {
    if (round(subclass) <= 1) {
      stop("subclass must be greater than 1.",call.=FALSE)
    }
  }
  else if (!all(subclass <= 1 & subclass >= 0)) {
    stop("When specifying subclass as a vector of quantiles, all values must be between 0 and 1.",
         call. = FALSE)
  }

  if (!is.null(sub.by)) {
    sub.by.choices <- c("treat", "control", "all")
    if (!is.character(sub.by) || length(sub.by) != 1 || anyNA(pmatch(sub.by, sub.by.choices))) {
      stop("'sub.by' is deprecated and can't be converted into a proper input. Please supply an argument to 'estimand' instead.", call. = FALSE)
    }
    else {
      sub.by <- sub.by.choices[pmatch(sub.by, sub.by.choices)]
      estimand <- switch(sub.by, "treat" = "ATT", "control" = "ATC", "ATE")
      warning(sprintf("'sub.by' is deprecated and has been replaced with 'estimand'. Setting 'estimand' to \"%s\".",
                      estimand), call. = FALSE, immediate. = TRUE)
    }
  }
  else {
    estimand <- toupper(estimand)
    estimand <- match_arg(estimand, c("ATT", "ATC", "ATE"))
  }

  if (is.null(min.n)) min.n <- 1
  else if (!is.numeric(min.n) || length(min.n) != 1) {
    stop("'min.n' must be a single number.", call. = FALSE)
  }

  n.obs <- length(treat)

  ## Setting Cut Points
  if (length(subclass) == 1) {
    sprobs <- seq(0, 1, length.out = round(subclass) + 1)
  }
  else {
    sprobs <- sort(subclass)
    if (sprobs[1] != 0) sprobs <- c(0, sprobs)
    if (sprobs[length(sprobs)] != 1) sprobs <- c(sprobs, 1)
    subclass <- length(sprobs) - 1
  }

  q <- switch(estimand,
              "ATT" = quantile(distance[treat==1], probs = sprobs, na.rm = TRUE),
              "ATC" = quantile(distance[treat==0], probs = sprobs, na.rm = TRUE),
              "ATE" = quantile(distance, probs = sprobs, na.rm = TRUE))

  ## Calculating Subclasses
  psclass <- setNames(rep(NA_integer_, n.obs), names(treat))
  psclass[!discarded] <- as.integer(findInterval(distance[!discarded], q, all.inside = TRUE))

  if (length(unique(na.omit(psclass))) != subclass){
    warning("Due to discreteness in the distance measure, fewer subclasses were generated than were requested.", call.=FALSE)
  }

  if (min.n == 0) {
    ## If any subclass are missing treated or control units, set all to NA
    is.na(psclass)[!discarded & !psclass %in% intersect(psclass[!discarded & treat == 1],
                                                        psclass[!discarded & treat == 0])] <- TRUE
  }
  else if (any(table(treat, psclass) < min.n)) {
    ## If any subclasses don't have members of a treatment group, fill them
    ## by "scooting" units from nearby subclasses until each subclass has a unit
    ## from each treatment group
    psclass[!discarded] <- subclass_scoot(psclass[!discarded], treat[!discarded], distance[!discarded], min.n)
  }

  psclass <- setNames(factor(psclass, nmax = length(q)), names(treat))
  levels(psclass) <- as.character(seq_len(nlevels(psclass)))

  if (verbose) cat("Calculating matching weights... ")

  res <- list(subclass = psclass, q.cut = q,
              weights = weights.subclass(psclass, treat, estimand))

  if (verbose) cat("Done.\n")

  class(res) <- c("matchit.subclass", "matchit")
  return(res)
}


# MATCHIT method = cardinality----------------------------------
matchit2cardinality <-  function(treat, data, discarded, formula,
                                 ratio = 1, focal = NULL, s.weights = NULL, replace = FALSE, exact = NULL,
                                 estimand = "ATT", verbose = FALSE,
                                 tols = .05, std.tols = TRUE,
                                 solver = "glpk", time = 1*60, ...){

  if (verbose) {
    cat("Cardinality matching... \n")
  }

  tvals <- unique(treat)
  nt <- length(tvals)

  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATT", "ATC", "ATE"))
  if (!is.null(focal)) {
    if (!focal %in% tvals) stop("'focal' must be a value of the treatment.", call. = FALSE)
  }
  else if (estimand == "ATC") {
    focal <- min(tvals)
  }
  else {
    focal <- max(tvals)
  }

  lab <- names(treat)

  weights <- setNames(rep(0, length(treat)), lab)

  X <- get.covs.matrix(formula, data = data)

  if (!is.null(exact)) {
    ex <- factor(exactify(model.frame(exact, data = data), nam = lab, sep = ", ", include_vars = TRUE))

    cc <- do.call("intersect", lapply(tvals, function(t) ex[treat == t]))
    if (length(cc) == 0) stop("No matches were found.", call. = FALSE)
  }
  else {
    ex <- NULL
  }

  #Process tols
  assign <- get_assign(X)

  if (length(tols) == 0 || !is.numeric(tols) || anyNA(tols)) stop("'tols' must be numeric.", call. = FALSE)
  if (length(tols) == 1) tols <- rep(tols, ncol(X))
  else if (length(tols) == max(assign)) {
    tols <- tols[assign]
  }
  else if (length(tols) != ncol(X)) {
    stop("'tols' must have length 1 or the number of covariates. See ?method_cardinality for details.", call. = FALSE)
  }

  if (length(std.tols) == 0 || !is.logical(std.tols) || anyNA(std.tols)) stop("'std.tols' must be logical (TRUE/FALSE).", call. = FALSE)
  if (length(std.tols) == 1) std.tols <- rep(std.tols, ncol(X))
  else if (length(std.tols) == max(assign)) {
    std.tols <- std.tols[assign]
  }
  else if (length(std.tols) != ncol(X)) {
    stop("'std.tols' must have length 1 or the number of covariates. See ?method_cardinality for details.", call. = FALSE)
  }

  #Apply std.tols
  if (any(std.tols)) {
    if (estimand == "ATE") {
      sds <- sqrt(Reduce("+", lapply(tvals, function(t) {
        apply(X[treat==t, std.tols, drop = FALSE], 2, wvar, w = s.weights[treat==t])
      }))/nt)
    }
    else {
      sds <- sqrt(apply(X[treat==focal,std.tols,drop=FALSE], 2, wvar, w = s.weights[treat==focal]))
    }

    zero.sds <- sds < 1e-10

    X[,std.tols][,!zero.sds] <- scale(X[, std.tols, drop = FALSE][,!zero.sds, drop = FALSE],
                                     center = FALSE, scale = sds[!zero.sds])
  }

  opt.out <- setNames(vector("list", if (is.null(ex)) 1L else nlevels(ex)), levels(ex))

  for (i in seq_along(opt.out)) {
    if (is.null(ex)) in.exact <- which(!discarded)
    else in.exact <- which(!discarded & ex == levels(ex)[i])

    out <- cardinality_matchit(treat = treat[in.exact],
                               X = X[in.exact,, drop = FALSE],
                               estimand = estimand, tols = tols,
                               s.weights = s.weights[in.exact],
                               ratio = ratio,
                               focal = focal, tvals = tvals,
                               solver = solver, time = time,
                               verbose = verbose)

    weights[in.exact] <- out[["weights"]]
    opt.out[[i]] <- out[["opt.out"]]
  }

  if (length(opt.out) == 1L) out <- out[[1]]

  psclass <- NULL

  res <- list(subclass = psclass,
              weights = weights,
              obj = opt.out)

  if (verbose) cat("Done.\n")

  class(res) <- "matchit"

  return(res)

}