# this function takes inputs from matchit() and returns the
# strata for each observation in the subclass entry, the
# weight for each observation in the weight entry, and the
# match.matrix object
#
# MATCHIT method= cem--------------------------------------
matchit2cem <- function(treat, covs, estimand = "ATT", verbose = FALSE, ...) {

  check.package("cem")

  if (verbose)
    cat("Coarsened exact matching...\n")

  A <- list(...)
  A[["method"]] <- A[["k2k.method"]]

  n.obs <- length(treat)

  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATT", "ATC", "ATE"))

  # cem takes the data all together and wants the treatment specified
  # with the column name of the data frame. Here we massage the matchit
  # inputs to this format. Note that X has its proper columnames, but
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
      stop(paste0("(from cem) ", conditionMessage(e)), call. = FALSE, immediate. = TRUE)
    }
  })

  # here we create a column vector where the matched entry get its stratum
  # and the unmatched entry gets an NA.
  strat <- setNames(rep(NA_character_, n.obs), names(treat))
  if (!is.null(mat)) strat[mat$matched] <- mat$strata[mat$matched]

  res <- list(subclass = strat,
              weights = weights.subclass(strat, treat, estimand))
  class(res) <- "matchit"
  return(res)
}

# MATCHIT method= exact------------------------------------
matchit2exact <- function(treat, covs, data, estimand = "ATT", verbose = FALSE, ...){

  if(verbose)
    cat("Exact matching... \n")

  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATT", "ATC", "ATE"))

  xx <- exactify(covs, names(treat))
  cc <- intersect(xx[treat==1], xx[treat==0])

  psclass <- setNames(match(xx, cc), names(treat))

  res <- list(subclass = psclass,
              weights = weights.subclass(psclass, treat, estimand))
  class(res) <- "matchit"
  return(res)
}

# MATCHIT method= full-------------------------------------
matchit2full <- function(treat, formula, data, distance, discarded,
                         caliper = NULL, mahvars = NULL, exact = NULL,
                         estimand = "ATT", verbose = FALSE,
                         is.full.mahalanobis, ...) {

  check.package("optmatch")

  if(verbose)
    cat("Full matching... \n")

  A <- list(...)

  fm.args <- c("min.controls", "max.controls", "omit.fraction", "mean.controls", "tol")
  A[!names(A) %in% fm.args] <- NULL

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

  if (sum(!discarded & treat != focal) < sum(!discarded & treat == focal)) {
    warning("Fewer ", tc[2], " units than ", tc[1], " units; some ", tc[2], " units will be matched to multiple ", tc[1], " units.", immediate. = TRUE, call. = FALSE)
  }

  treat_ <- setNames(as.integer(treat == focal), names(treat))
  treat_[discarded] <- NA

  within.match <- NULL
  if (!is.null(exact)) {
    environment(exact) <- sys.frame(sys.nframe())
    within.match <- optmatch::exactMatch(update(exact, treat_ ~ .), data = data)
  }

  if (!is.null(caliper)) {
    if (any(names(caliper) != "")) {
      cov.cals <- setdiff(names(caliper), "")
      calcovs <- get.covs.matrix(reformulate(cov.cals, intercept = FALSE), data = data)
    }
    for (i in seq_along(caliper)) {
      if (names(caliper)[i] == "") {
        mo <- optmatch::match_on(distance, z = treat_)
      }
      else {
        mo <- optmatch::match_on(setNames(calcovs[,names(caliper)[i]], names(treat)), z = treat_)
      }

      if (is.null(within.match)) within.match <- optmatch::caliper(mo, caliper[i])
      else within.match <- within.match + optmatch::caliper(mo, caliper[i])

      rm(mo)
    }
  }

  withCallingHandlers({
    if (is.full.mahalanobis) {
      formula <- update(formula, treat_ ~ .)
      environment(formula) <- sys.frame(sys.nframe())

      full <- do.call(optmatch::fullmatch,
                      c(list(formula,
                             data = data,
                             method = "mahalanobis",
                             within = within.match),
                        A))
    }
    else if (!is.null(mahvars)) {
      mahvars <- update(mahvars, treat_ ~ .)
      environment(mahvars) <- sys.frame(sys.nframe())

      full <- do.call(optmatch::fullmatch,
                      c(list(mahvars,
                             data = data,
                             method = "mahalanobis",
                             within = within.match),
                        A))
    }
    else {
      full <- do.call(optmatch::fullmatch,
                      c(list(treat_ ~ distance,
                             method = "euclidean", #slightly faster than Mahalanobis
                             within = within.match),
                        A))
    }
  },
  warning = function(w) {
    warning(paste0("(from optmatch) ", conditionMessage(w)), call. = FALSE, immediate. = TRUE)
    invokeRestart("muffleWarning")
  },
  error = function(e) {
    stop(paste0("(from optmatch) ", conditionMessage(e)), call. = FALSE, immediate. = TRUE)
  })

  if (all(is.na(full))) stop("No matches were found.", call. = FALSE)

  psclass <- as.integer(as.factor(full))
  names(psclass) <- names(treat)

  #No match.matrix because treated units don't index matched strata (i.e., more than one
  #treated unit can be in the same stratum). Stratum information is contained in subclass.

  res <- list(subclass = psclass,
              weights = weights.subclass(psclass, treat, estimand))

  class(res) <- c("matchit")
  return(res)
}

# MATCHIT method= optimal----------------------------------
matchit2optimal <- function(treat, formula, data, distance, discarded,
                            ratio = 1, caliper = NULL, mahvars = NULL, exact = NULL,
                            estimand = "ATT", verbose = FALSE,
                            is.full.mahalanobis,...) {

  check.package("optmatch")

  if (verbose) cat("Optimal matching... \n")

  A <- list(...)
  pm.args <- c("tol")
  A[!names(A) %in% pm.args] <- NULL

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

  if (sum(!discarded & treat != focal) < sum(!discarded & treat == focal)) {
    warning("Fewer ", tc[2], " units than ", tc[1], " units; not all ", tc[1], " units will get a match.", immediate. = TRUE, call. = FALSE)    }
  else if (sum(!discarded & treat != focal) < sum(!discarded & treat == focal)*ratio) {
    stop(paste0("Not enough ", tc[2], " units for ", ratio, " matches for each ", tc[1], " unit."), call. = FALSE)
  }

  treat_ <- setNames(as.integer(treat == focal), names(treat))
  treat_[discarded] <- NA

  ratio <- process.ratio(ratio)

  if (!is.null(exact)) {
    environment(exact) <- sys.frame(sys.nframe())
    exact.match <- optmatch::exactMatch(update(exact, treat_ ~ .), data = data)
  }
  else exact.match <- NULL

  if (!is.null(caliper)) {
    warning("Calipers are currently not compatible with method = \"optimal\" and will be ignored.", call. = FALSE)
    caliper <- NULL
  }

  withCallingHandlers({
    if (is.full.mahalanobis) {
      formula <- update(formula, treat_ ~ .)
      environment(formula) <- sys.frame(sys.nframe())

      pair <- do.call(optmatch::pairmatch,
                      c(list(formula,
                             data = data,
                             method = "mahalanobis",
                             controls = ratio,
                             within = exact.match),
                        A))
    }
    else if (!is.null(mahvars)) {
      mahvars <- update(mahvars, treat_ ~ .)
      environment(mahvars) <- sys.frame(sys.nframe())

      pair <- do.call(optmatch::pairmatch,
                      c(list(mahvars,
                             data = data,
                             method = "mahalanobis",
                             controls = ratio,
                             within = exact.match),
                        A))
    }
    else {
      pair <- do.call(optmatch::pairmatch,
                      c(list(treat_ ~ distance,
                             method = "euclidean", #slightly faster than Mahalanobis
                             controls = ratio,
                             within = exact.match),
                        A))
    }
  },
  warning = function(w) {
    warning(paste0("(from optmatch) ", conditionMessage(w)), call. = FALSE, immediate. = TRUE)
    invokeRestart("muffleWarning")
  },
  error = function(e) {
    stop(paste0("(from optmatch) ", conditionMessage(e)), call. = FALSE, immediate. = TRUE)
  })

  if (all(is.na(pair))) stop("No matches were found.", call. = FALSE)

  psclass <- as.integer(as.factor(pair))
  names(psclass) <- names(treat)
  na.class <- is.na(psclass)

  mm <- matrix(NA_character_, ncol = ratio, nrow = sum(treat == focal),
               dimnames = list(names(treat)[treat == focal], seq_len(ratio)))

  for (i in rownames(mm)[!na.class[treat == focal]]) {
    matched.units <- names(treat)[treat != focal & !na.class & psclass == psclass[i]]
    if (length(matched.units) > 0) mm[i, seq_along(matched.units)] <- matched.units
  }

  ## calculate weights and return the results
  res <- list(match.matrix = mm, subclass = psclass,
              weights = weights.subclass(psclass, treat, estimand))

  class(res) <- "matchit"
  return(res)
}
# MATCHIT method= genetic----------------------------------
#Needs updates
matchit2genetic <- function(treat, data, distance, discarded,
                            ratio = 1, replace = FALSE, m.order = NULL,
                            caliper = NULL, mahvars = NULL, exact = NULL,
                            formula = NULL, estimand = "ATT", verbose = FALSE,
                            is.full.mahalanobis, use.genetic = TRUE, ...) {

  check.package(c("Matching", "rgenoud"))

  if (verbose)
    cat("Genetic matching... \n")

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
      warning("Fewer ", tc[2], " units than ", tc[1], " units; not all ", tc[1], " units will get a match.", immediate. = TRUE, call. = FALSE)    }
    else if (sum(!discarded & treat != focal) < sum(!discarded & treat == focal)*ratio) {
      stop(paste0("Not enough ", tc[2], " units for ", ratio, " matches for each ", tc[1], " unit."), call. = FALSE)
    }
  }

  treat <- setNames(as.integer(treat == focal), names(treat))

  ratio <- process.ratio(ratio)

  n.obs <- length(treat)
  n1 <- sum(treat == 1)
  n0 <- sum(treat == 0)

  if (is.null(names(treat))) names(treat) <- seq_len(n.obs)

  if (!is.null(distance)) {
    if (is.null(m.order)) m.order <- if (estimand == "ATC") "smallest" else "largest"
    else m.order <- match_arg(m.order, c("largest", "smallest", "random", "data"))
    ord <- switch(m.order,
                  "largest" = order(distance, decreasing = TRUE),
                  "smallest" = order(distance),
                  "random" = sample(seq_len(n.obs), n.obs, replace = FALSE),
                  "data" = seq_len(n.obs))
  }
  else {
    m.order <- match_arg(m.order, c("data", "random"))
    ord <- switch(m.order,
                  "random" = sample(seq_len(n.obs), n.obs, replace = FALSE),
                  "data" = seq_len(n.obs))
  }
  ord <- ord[!ord %in% which(discarded)]

  #Create X (matching variables) and covs_to_balance
  if (!is.null(mahvars)) {
    covs_to_balance <- get.covs.matrix(formula, data = data)
    X <- get.covs.matrix(mahvars, data = data)
  }
  else if (is.full.mahalanobis) {
    covs_to_balance <- get.covs.matrix(formula, data = data)
    X <- covs_to_balance
  }
  else {
    covs_to_balance <- get.covs.matrix(formula, data = data)
    X <- cbind(covs_to_balance, distance)
  }

  #Process exact; exact.log will be supplied to GenMatch() and Match()
  if (!is.null(exact)) {
    #Add covariates in exact not in X to X
    ex <- as.numeric(factor(exactify(model.frame(exact, data = data), names(treat))))
    X <- cbind(X, ex)

    exact.log <- c(rep(FALSE, ncol(X) - 1), TRUE)
  }
  else exact.log <- NULL

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
      if (!is.null(mahvars)) {
        #If mahvars specified, distance is not yet in X, so add it to X
        X <- cbind(X, distance)
        cal <- c(cal, caliper[names(caliper) == ""])

        #Expand exact.log for newly added distance
        if (!is.null(exact.log)) exact.log <- c(exact.log, FALSE)
      }
      else {
        #Otherwise, distance is in X at the specified index
        cal[ncol(covs_to_balance) + 1] <- caliper[names(caliper) == ""]
      }
    }

  }
  else cal<- NULL

  #Reorder data according to m.order since Match matches in order of data;
  #ord already excludes discarded units

  treat_ <- treat[ord]
  covs_to_balance <- covs_to_balance[ord,,drop = FALSE]
  X <- X[ord,,drop = FALSE]

  if (use.genetic) {
    withCallingHandlers({
      g.out <- do.call(Matching::GenMatch,
                       c(list(Tr = treat_, X = X, BalanceMatrix = covs_to_balance,
                              M = ratio, exact = exact.log, caliper = cal,
                              replace = replace, estimand = "ATT", ties = FALSE,
                              CommonSupport = FALSE, verbose = verbose,
                              print.level = 2*verbose),
                         A[names(A) %in% names(formals(Matching::GenMatch))]))
    },
    warning = function(w) {
      if (!startsWith(conditionMessage(w), "replace==FALSE, but there are more (weighted) treated obs than control obs.")) {
        warning(paste0("(from Matching) ", conditionMessage(w)), call. = FALSE, immediate. = TRUE)
      }
      invokeRestart("muffleWarning")
    },
    error = function(e) {
      stop(paste0("(from Matching) ", conditionMessage(e)), call. = FALSE, immediate. = TRUE)
    })
  }
  else {
    #For debugging
    g.out <- NULL
  }

  withCallingHandlers({
    m.out <- Matching::Match(Tr = treat_, X = X,
                             M = ratio, exact = exact.log, caliper = cal,
                             replace = replace, estimand = "ATT", ties = FALSE,
                             CommonSupport = FALSE, Weight = if (use.genetic) 3 else 2,
                             Weight.matrix = g.out, version = "fast")
  },
  warning = function(w) {
    if (!startsWith(conditionMessage(w), "replace==FALSE, but there are more (weighted) treated obs than control obs.")) {
      warning(paste0("(from Matching) ", conditionMessage(w)), call. = FALSE, immediate. = TRUE)
    }
    invokeRestart("muffleWarning")
  },
  error = function(e) {
    stop(paste0("(from Matching) ", conditionMessage(e)), call. = FALSE, immediate. = TRUE)
  })

  lab <- names(treat)
  lab1 <- names(treat[treat == 1])

  lab_ <- names(treat_)

  mm <- matrix(NA_character_, nrow = n1, ncol = max(table(m.out$index.treated)),
               dimnames = list(lab1, NULL))

  unique.matched.focal <- unique(m.out$index.treated, nmax = n1)

  for (i in unique.matched.focal) {
    matched.units <- lab_[m.out$index.control[m.out$index.treated == i]]
    mm[lab_[i], seq_along(matched.units)] <- matched.units
  }

  if (!replace) {
    psclass <- setNames(rep(NA_character_, n.obs), lab)
    no.match <- is.na(mm)
    psclass[lab1[!no.match[,1]]] <- lab1[!no.match[,1]]
    psclass[mm[!no.match]] <- lab1[row(mm)[!no.match]]
    psclass <- setNames(as.integer(factor(psclass, nmax = n1)), lab)
  }
  else psclass <- NULL

  res <- list(match.matrix = mm,
              subclass = psclass,
              weights = weights.matrix(mm, treat),
              obj = g.out)

  class(res) <- "matchit"
  return(res)
}

# MATCHIT method= nearest----------------------------------
matchit2nearest <-  function(treat, data, distance, discarded,
                             ratio = 1, replace = FALSE, m.order = NULL,
                             caliper = NULL, mahvars = NULL, exact = NULL,
                             formula = NULL, estimand = "ATT", verbose = FALSE,
                             is.full.mahalanobis, ...){

  if(verbose)
    cat("Nearest neighbor matching... \n")

  ratio <- process.ratio(ratio)

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
      warning("Fewer ", tc[2], " units than ", tc[1], " units; not all ", tc[1],
              " units will get a match.", immediate. = TRUE, call. = FALSE)    }
    else if (sum(!discarded & treat != focal) < sum(!discarded & treat == focal)*ratio) {
      warning(paste0("Not enough ", tc[2], " units for ", ratio, " matches for each ", tc[1],
                     " unit. Some ", tc[1], " units will not be matched to ", ratio, " ", tc[2], " units."), immediate. = TRUE, call. = FALSE)
    }
  }

  treat <- setNames(as.integer(treat == focal), names(treat))

  n.obs <- length(treat)
  n1 <- sum(treat == 1)
  n0 <- sum(treat == 0)

  dis1 <- discarded[treat == 1]
  dis0 <- discarded[treat == 0]

  n1_ <- sum(!dis1)
  n0_ <- sum(!dis0)

  lab1 <- names(treat[treat == 1])
  lab0 <- names(treat[treat == 0])

  if (!is.null(distance)) {
    names(distance) <- names(treat)
    d1 <- distance[treat == 1]
    d0 <- distance[treat == 0]
  }

  if (!is.null(caliper)) {
    if (any(names(caliper) != "")) {
      calcovs <- get.covs.matrix(reformulate(setdiff(names(caliper), "")), data = data)
      rownames(calcovs) <- names(treat)
    }

    caliper.dist <- caliper[names(caliper) == ""]
    caliper <- caliper[names(caliper) != ""]
  }
  else {
    caliper.dist <- NULL
  }

  if (!is.null(exact)) {
    ex <- exactify(model.frame(exact, data = data), names(treat))
    ex1 <- ex[treat == 1]
    ex0 <- ex[treat == 0]
  }

  if (is.full.mahalanobis) {
    mahcovs <- get.covs.matrix(formula, data)
    mahSigma_inv <- generalized_inverse(cov(mahcovs))
    rownames(mahcovs) <- names(treat)
  }
  else if (!is.null(mahvars)) {
    mahcovs <- get.covs.matrix(mahvars, data)
    mahSigma_inv <- generalized_inverse(cov(mahcovs))
    rownames(mahcovs) <- names(treat)
  }
  else {
    mahcovs <- mahSigma_inv <- NULL
  }

  mm <- matrix(NA_character_, nrow = n1,
               ncol = ratio, dimnames = list(lab1, seq_len(ratio)))

  unmatched <- setNames(rep(TRUE, n0_), lab0[!dis0])

  if (replace) {
    ord <- seq_len(n1)
  }
  else {
    if (!is.null(distance)) {
      if (is.null(m.order)) m.order <- if (estimand == "ATC") "smallest" else "largest"
      else m.order <- match_arg(m.order, c("largest", "smallest", "random", "data"))
      ord <- switch(m.order,
                    "largest" = order(d1, decreasing = TRUE),
                    "smallest" = order(d1),
                    "random" = sample(seq_len(n1), n1, replace = FALSE),
                    "data" = seq_len(n1))
    }
    else {
      m.order <- match_arg(m.order, c("data", "random"))
      ord <- switch(m.order,
                    "random" = sample(seq_len(n1), n1, replace = FALSE),
                    "data" = seq_len(n1))
    }
  }

  for (r in seq_len(ratio)) {
    for (i in seq_len(n1)[!dis1]) {

      ord_i <- ord[i]
      c.eligible <- names(unmatched)

      if (!replace) {
        if (!any(unmatched)) break
        c.eligible <- c.eligible[unmatched]
      }
      else if (r > 1) {
        #If replace = T and r > 1, don't rematch to same control unit
        c.eligible <- c.eligible[!c.eligible %in% mm[ord_i, seq_len(r - 1)]]
      }

      if (!is.null(exact) && length(c.eligible) > 0) {
        c.eligible <- c.eligible[ex0[c.eligible] == ex1[ord_i]]
      }

      #Get distances among eligible and apply caliper
      if (length(c.eligible) > 0) {

        ps.diff <- NULL

        #PS caliper
        if (length(caliper.dist) > 0) {
          ps.diff <- abs(unname(d1[ord_i]) - d0[c.eligible])
          c.eligible <- c.eligible[ps.diff <= caliper.dist]
        }

        #Covariate caliper
        if (length(caliper) > 0 && length(c.eligible) > 0) {
          for (x in names(caliper)) {
            calcov.diff <- abs(unname(calcovs[lab1[ord_i], x]) - calcovs[c.eligible, x])
            c.eligible <- c.eligible[calcov.diff <= caliper[x]]
            if (length(c.eligible) == 0) break
          }
        }

        if (is.null(mahcovs)) {
          #PS matching
          distances <- if(is.null(ps.diff)) abs(unname(d1[ord_i]) - d0[c.eligible]) else ps.diff
        }
        else {
          #MD matching
          distances <- sqrt(mahalanobis(mahcovs[c.eligible, ,drop = FALSE],
                                        mahcovs[lab1[ord_i],],
                                        cov = mahSigma_inv, inverted = TRUE))
        }
        # names(distances) <- c.eligible
      }

      #Assign match
      if (length(c.eligible) > 0) {

        min.dist <- min(distances[c.eligible])

        #For now, resolve ties by selecting the first unit
        mm[ord_i, r] <- c.eligible[distances[c.eligible] == min.dist][1]

        if (!replace) unmatched[mm[ord_i, r]] <- FALSE

      }
    }
  }

  if (!replace) {
    psclass <- setNames(rep(NA_character_, n.obs), names(treat))
    no.match <- is.na(mm)
    psclass[lab1[!no.match[,1]]] <- lab1[!no.match[,1]]
    psclass[mm[!no.match]] <- lab1[row(mm)[!no.match]]
    psclass <- setNames(as.integer(factor(psclass, nmax = n1)), names(treat))
  }
  else psclass <- NULL

  res <- list(match.matrix = mm,
              subclass = psclass,
              weights = weights.matrix(mm, treat))

  class(res) <- "matchit"

  return(res)
}

# MATCHIT method= subclass---------------------------------
#Needs updates
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
  else if (!is.vector(subclass, "numeric")) {
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
    if (!is.vector(sub.by, "character") || length(sub.by) != 1 || anyNA(pmatch(sub.by, sub.by.choices))) {
      stop("'sub.by' is deprecated and can't be converted into a proper input. Please supply an argument to 'estimand' instead.", call. = FALSE)
    }
    else {
      sub.by <- sub.by.choices[pmatch(sub.by, sub.by.choices)]
      estimand <- switch(sub.by, "treat" = "ATT", "control" = "ATC", "ATE")
      warning(paste0("'sub.by' is deprecated and has been replaced with 'estimand'. Setting 'estimand' to \"", estimand, "\"."), call. = FALSE, immediate. = TRUE)
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

  p1 <- distance[treat==1]
  p0 <- distance[treat==0]

  ## Setting Cut Points
  if (length(subclass) == 1) {
    sprobs <- seq(0, 1, length = round(subclass) + 1)
  }
  else {
    sprobs <- sort(subclass)
    if (sprobs[1] != 0) sprobs <- c(0, sprobs)
    if (sprobs[length(sprobs)] != 1) sprobs <- c(sprobs, 1)
    subclass <- length(sprobs) - 1
  }

  q <- switch(estimand,
              "ATT" = quantile(p1, probs = sprobs, na.rm = TRUE),
              "ATC" = quantile(p0, probs = sprobs, na.rm = TRUE),
              "ATE" = quantile(distance, probs = sprobs, na.rm = TRUE))

  ## Calculating Subclasses
  psclass <- setNames(rep(NA_integer_, n.obs), names(treat))
  psclass[!discarded] <- as.integer(findInterval(distance[!discarded], q, all.inside = TRUE))

  ## If any subclasses don't have members of a treatment group, fill them
  ## by "scooting" units from nearby subclasses until each subclass has a unit
  ## from each treatment group
  if (any(table(treat, psclass) < min.n)) {
    psclass[!discarded] <- subclass_scoot(psclass[!discarded], treat[!discarded], distance[!discarded], min.n)
  }

  if (verbose) cat("Done\n")

  #warning for discrete data
  unique.classes <- unique(psclass[!is.na(psclass)], nmax = subclass)

  if (length(unique.classes) != subclass){
    warning("Due to discreteness in the distance measure, fewer subclasses were generated than were requested.", call.=FALSE)
  }

  res <- list(subclass = psclass, q.cut = q,
              weights = weights.subclass(psclass, treat, estimand))

  class(res) <- c("matchit.subclass", "matchit")
  return(res)
}

