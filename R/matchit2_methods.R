# this function takes inputs from matchit() and returns the
# strata for each observation in the subclass entry and the
# weight for each observation in the weight entry. No match
# matrix is returned since matches are not unique within
# strata.
#
# MATCHIT method= cem--------------------------------------
matchit2cem <- function(treat, covs, data, distance, discarded, is.full.mahalanobis,
                            ratio = 1, verbose = FALSE, ...) {

  check.package("cem")

  if (verbose)
    cat("Coarsened exact matching...\n")

  A <- list(...)
  A[["method"]] <- A[["k2k.method"]]

  n.obs <- length(treat)

  # cem takes the data all together and wants the treatment specified
  # with the column name of the data frame. Here we massage the matchit
  # inputs to this format. Note that X has its proper columnames, but
  # treat does not have the original column name.
  cem.data <- data.frame(covs, treat)

  args.excluded <- c("treatment", "baseline.group", "data", "verbose")
  mat <- do.call(cem::cem, c(list(treatment = names(cem.data)[ncol(cem.data)],
                                  data = cem.data,
                                  verbose = as.integer(verbose)),
                             A[names(A) %in% setdiff(names(formals(cem::cem)), args.excluded)]))

  # here we create a column vector where the matched entry get its stratum
  # and the unmatched entry gets an NA.
  strat <- setNames(rep(NA_character_, n.obs), names(treat))
  strat[mat$matched] <- mat$strata[mat$matched]

  # here we just add the names onto the wieght from the cem output
  wh <- setNames(mat$w, names(treat))

  # weighting functions in matchit error-out on these conditions,
  # so we should too.

  if (sum(wh)==0)
    stop("No units were matched")
  else if (sum(wh[treat==1])==0)
    stop("No treated units were matched")
  else if (sum(wh[treat==0])==0)
    stop("No control units were matched")

  res <- list(subclass = strat, weights = mat$w)
  class(res) <- "matchit"
  return(res)
}

# MATCHIT method= exact------------------------------------
matchit2exact <- function(treat, covs, data, distance, discarded, is.full.mahalanobis, verbose=FALSE, ...){

  if(verbose)
    cat("Exact matching... \n")

  xx <- setNames(do.call("paste", c(covs, sep = "|")), names(treat))
  cc <- intersect(xx[treat==1], xx[treat==0])

  psclass <- setNames(match(xx, cc), names(treat))

  res <- list(subclass = psclass, weights = weights.subclass(psclass, treat))
  class(res) <- c("matchit.exact", "matchit")
  return(res)
}

# MATCHIT method= full-------------------------------------
matchit2full <- function(treat, covs, data, distance, discarded,
                         ratio = 1, caliper = NULL, mahvars = NULL, exact = NULL,
                         formula = NULL, verbose = FALSE,
                         is.full.mahalanobis, ...) {

  check.package("optmatch")

  if(verbose)
    cat("Full matching... \n")

  treat0 <- treat
  treat0[discarded] <- NA

  if (sum(!discarded & treat == 0) < sum(!discarded & treat == 1)) {
    warning("Fewer control units than treated units; some control units will be matched to multiple treated units.", immediate. = TRUE, call. = FALSE)
  }

  if (!is.null(exact)) {
    ex <- do.call("paste", c(model.frame(exact, data = data), sep = "|"))
    exact.match <- optmatch::exactMatch(ex, treat0)
  }
  else exact.match <- NULL

  if (is.full.mahalanobis) {
    data <- data.frame(treat0, covs)
    full <- optmatch::fullmatch(formula(data),
                                data = data,
                                method = "mahalanobis",
                                within = exact.match,
                                ...)
  }
  else if (!is.null(mahvars)) {
    full <- optmatch::fullmatch(update(mahvars, treat0 ~ .),
                                data = data,
                                method = "mahalanobis",
                                within = exact.match,
                                ...)
  }
  else {
    full <- optmatch::fullmatch(treat0 ~ distance,
                                method = "euclidean", #slightly faster than Mahalanobis
                                within = exact.match,
                                ...)
  }

  psclass <- as.integer(as.factor(full))
  names(psclass) <- names(treat)
  na.class <- is.na(psclass)
  max.matches <- max(vapply(levels(full), function(x) sum(!na.class & full == x) - 1, numeric(1L)))

  mm <- matrix(NA_character_, ncol = max.matches, nrow = sum(treat == 1),
               dimnames = list(names(treat)[treat == 1], 1:max.matches))

  for (i in rownames(mm)[!is.na(psclass[treat == 1])]) {
    matched.controls <- names(treat)[treat == 0 & !na.class & psclass == psclass[i]]
    if (length(matched.controls) > 0) mm[i, 1:length(matched.controls)] <- matched.controls
  }

  res <- list(match.matrix = mm, subclass = psclass,
              weights = weights.subclass(psclass, treat)
              # weights = weights.matrix(mm, treat, discarded) #slower
              )

  if (is.full.mahalanobis) class(res) <- c("matchit.mahalanobis", "matchit")
  else class(res) <- c("matchit")
  return(res)
}

# MATCHIT method= optimal----------------------------------
matchit2optimal <- function(treat, covs, data, distance, discarded,
                            ratio = 1, caliper = NULL, mahvars = NULL, exact = NULL,
                            formula = NULL, verbose = FALSE,
                            is.full.mahalanobis,...) {

  check.package("optmatch")

  if(verbose)
    cat("Optimal matching... \n")

  treat0 <- treat
  treat0[discarded] <- NA

  ratio <- process.ratio(ratio)

  if (sum(!discarded & treat == 0) < sum(!discarded & treat == 1)) {
    warning("Fewer control units than treated units; not all treated units will get a match.", immediate. = TRUE, call. = FALSE)
  }
  else if (sum(!discarded & treat == 0) < sum(!discarded & treat == 1)*ratio) {
    stop(paste0("Not enough control units for ", ratio, " matches for each treated unit."), call. = FALSE)
  }

  if (!is.null(exact)) {
    exact.match <- optmatch::exactMatch(update(exact, treat0 ~ .), data = data)
  }
  else exact.match <- NULL

  if (!is.null(caliper)) {
    warning("Calipers are currently not compatible with method = \"optimal\" and will be ignored.", call. = FALSE)
    caliper <- NULL
  }

  if (is.full.mahalanobis) {
    data <- data.frame(treat0, covs)
    pair <- optmatch::pairmatch(formula(data),
                                data = data,
                                method = "mahalanobis",
                                controls = ratio,
                                within = exact.match,
                                ...)
  }
  else if (!is.null(mahvars)) {
    pair <- optmatch::pairmatch(update(mahvars, treat0 ~ .),
                                data = data,
                                method = "mahalanobis",
                                controls = ratio,
                                within = exact.match,
                                ...)
  }
  else {
    pair <- optmatch::pairmatch(treat0 ~ distance,
                                method = "euclidean", #slightly faster than Mahalanobis
                                controls = ratio,
                                within = exact.match, ...)
  }

  psclass <- as.integer(as.factor(pair))
  names(psclass) <- names(treat)
  na.class <- is.na(psclass)

  mm <- matrix(NA_character_, ncol = ratio, nrow = sum(treat == 1),
               dimnames = list(names(treat)[treat == 1], seq_len(ratio)))

  for (i in rownames(mm)[!na.class[treat == 1]]) {
    matched.controls <- names(treat)[treat == 0 & !na.class & psclass == psclass[i]]
    if (length(matched.controls) > 0) mm[i, seq_along(matched.controls)] <- matched.controls
  }

  ## calculate weights and return the results
  res <- list(match.matrix = mm, subclass = psclass,
              weights = weights.matrix(mm, treat, discarded))

  if (is.full.mahalanobis) class(res) <- c("matchit.mahalanobis", "matchit")
  else class(res) <- "matchit"
  return(res)
}
# MATCHIT method= genetic----------------------------------
#Needs updates
matchit2genetic <- function(treat, covs, data, distance, discarded,
                            ratio = 1, replace = FALSE, m.order = NULL,
                            caliper = NULL, mahvars = NULL, exact = NULL,
                            formula = NULL, verbose = FALSE,
                            is.full.mahalanobis, use.genetic = TRUE, ...) {

  check.package(c("Matching", "rgenoud"))

  if (verbose)
    cat("Genetic matching... \n")

  A <- list(...)

  ratio <- process.ratio(ratio)

  n.obs <- length(treat)
  n1 <- sum(treat == 1)
  n0 <- sum(treat == 0)

  if (!replace) {
    if (sum(!discarded & treat == 0) < sum(!discarded & treat == 1)) {
      warning("Fewer control units than treated units; not all treated units will get a match.", immediate. = TRUE, call. = FALSE)
    }
    else if (sum(!discarded & treat == 0) < sum(!discarded & treat == 1)*ratio) {
      warning(paste0("Not enough control units for ", ratio, " matches for each treated unit. Some treated units will not be matched."), immediate. = TRUE, call. = FALSE)
    }
  }

  if (is.null(names(treat))) names(treat) <- 1:n.obs

  if (!is.null(distance)) {
    m.order <- match_arg(m.order, c("largest", "smallest", "random", "data"))
    ord <- switch(m.order,
                  "largest" = order(distance, decreasing = TRUE),
                  "smallest" = order(distance),
                  "random" = sample(1:n.obs, n.obs, replace = FALSE),
                  "data" = 1:n.obs)
  }
  else {
    m.order <- match_arg(m.order, c("data", "random"))
    ord <- switch(m.order,
                  "random" = sample(1:n.obs, n.obs, replace = FALSE),
                  "data" = 1:n.obs)
  }

  if (!is.null(mahvars)) {
    covs_to_balance <- model.matrix(update(formula, NULL ~ . + 1), data = covs)[,-1, drop = FALSE]
    X <- cbind(distance, model.matrix(update(mahvars, .~.+1), data = data))[,-1, drop = FALSE]
  }
  else if (is.full.mahalanobis) {
    covs_to_balance <- model.matrix(update(formula, NULL ~ . + 1), data = covs)[,-1, drop = FALSE]
    X <- covs_to_balance
  }
  else {
    covs_to_balance <- model.matrix(update(formula, NULL ~ . + 1), data = covs)[,-1, drop = FALSE]
    X <- cbind(distance, covs_to_balance)
  }

  if (!is.null(exact)) {
    ex <- setNames(as.numeric(factor(do.call("paste", c(model.frame(exact, data = data), sep = "|")))), names(treat))
    X <- cbind(X, ex)
    exact.log <- c(rep(FALSE, ncol(X) - 1), TRUE)
  }
  else exact.log <- NULL

  if (!is.null(caliper)) {
    if (is.full.mahalanobis) {
      stop("A caliper cannot be set for genetic matching with distance = \"mahalanobis\".", call. = FALSE)
    }
    caliper <- c(caliper, rep(Inf, ncol(X) - 1))
  }
  else caliper <- NULL

  lab <- names(treat)
  lab1 <- names(treat[treat == 1])
  treat <- treat

  #Reorder data according to m.order since Match matches in order of data
  olab <- names(treat[ord])
  olab1 <- names(treat[ord][treat[ord] == 1])

  otreat <- treat[ord][!discarded[ord]]
  ocovs_to_balance <- covs_to_balance[ord,,drop = FALSE][!discarded[ord],, drop = FALSE]
  oX <- X[ord,,drop = FALSE][!discarded[ord],, drop = FALSE]

  if (use.genetic) {
    withCallingHandlers({
      g.out <- Matching::GenMatch(otreat, X = oX, BalanceMatrix = ocovs_to_balance,
                                  M = ratio, exact = exact.log, caliper = caliper,
                                  replace = replace, estimand = "ATT",
                                  CommonSupport = FALSE,
                                  verbose = verbose, print.level = 2*verbose, ...)
    },
    warning = function(w) {
      if (!startsWith(conditionMessage(w), "replace==FALSE, but there are more (weighted) treated obs than control obs.")) warning(w)
      invokeRestart("muffleWarning")
    })
  }
  else {
    #For debugging
    g.out <- NULL
  }

  withCallingHandlers({
    m.out <- Matching::Match(Tr = otreat, X = oX,
                             M = ratio, exact = exact.log, caliper = caliper,
                             replace = replace, estimand = "ATT",
                             CommonSupport = FALSE, Weight = if (use.genetic) 3 else 2,
                             Weight.matrix = g.out, version = "fast", ...)
  },
  warning = function(w) {
    if (!startsWith(conditionMessage(w), "replace==FALSE, but there are more (weighted) treated obs than control obs.")) warning(w)
    invokeRestart("muffleWarning")
  })

  mm <- matrix(NA_character_, nrow = n1, ncol = max(table(m.out$index.treated)),
               dimnames = list(lab1, NULL))

  unique.matched.treated <- unique(m.out$index.treated, nmax = n1)
  unique.matched.control <- unique(m.out$index.control, nmax = n0)

  for (i in unique.matched.treated) {
    matched.controls <- olab[m.out$index.control[m.out$index.treated == i]]
    mm[olab[i], 1:length(matched.controls)] <- matched.controls
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
              weights = weights.matrix(mm, treat, discarded))

  class(res) <- "matchit"
  return(res)
}

# MATCHIT method= nearest----------------------------------
matchit2nearest <-  function(treat, data, distance, discarded,
                             ratio = 1, replace = FALSE, m.order = NULL,
                             caliper = NULL, mahvars = NULL, exact = NULL,
                             formula = NULL, verbose = FALSE,
                             is.full.mahalanobis, ...){

  if(verbose)
    cat("Nearest neighbor matching... \n")

  ratio <- process.ratio(ratio)

  n.obs <- length(treat)
  n1 <- sum(treat == 1)
  n0 <- sum(treat == 0)

  if (!replace) {
    if (n0 < n1) {
      warning("Fewer control units than treated units; not all treated units will get a match.", immediate. = TRUE, call. = FALSE)
    }
    else if (n0 < n1*ratio) {
      warning(paste0("Not enough control units for ", ratio, " matches for each treated unit. Some treated units will receive fewer than ", ratio, " matches."), immediate. = TRUE, call. = FALSE)
    }
  }

  if (is.null(names(treat))) names(treat) <- 1:n.obs

  lab1 <- names(treat[treat == 1])
  lab0 <- names(treat[treat == 0])

  if (!is.null(caliper) || !is.full.mahalanobis) {
    if (is.null(names(distance))) names(distance) <- names(treat)
    cal <- caliper*sd(distance[!discarded])
    d1 <- distance[treat == 1]
    d0 <- distance[treat == 0]
  }

  dis1 <- discarded[treat == 1]
  dis0 <- discarded[treat == 0]

  if (!is.null(exact)) {
    ex <- setNames(do.call("paste", c(model.frame(exact, data = data), sep = "|")), names(treat))
    ex1 <- ex[treat == 1]
    ex0 <- ex[treat == 0]
  }

  if (is.full.mahalanobis) {
    mahcovs <- model.matrix(update(formula, NULL ~ . + 1), data)[,-1,drop = FALSE]
    mahSigma_inv <- MASS::ginv(cov(mahcovs))
  }
  else if (!is.null(mahvars)) {
    mahcovs <- model.matrix(update(mahvars, NULL ~ . + 1), data)[,-1,drop = FALSE]
    mahSigma_inv <- MASS::ginv(cov(mahcovs))
  }
  else mahcovs <- NULL

  mm <- matrix(NA_character_, nrow = n1,
                         ncol = ratio, dimnames = list(lab1, 1:ratio))

  unmatched <- setNames(rep(TRUE, sum(!dis0)), lab0[!dis0])

  if (!is.null(distance)) {
    m.order <- match_arg(m.order, c("largest", "smallest", "random", "data"))
    ord <- switch(m.order,
                  "largest" = order(d1, decreasing = TRUE),
                  "smallest" = order(d1),
                  "random" = sample(1:n1, n1, replace = FALSE),
                  "data" = 1:n1)
  }
  else {
    m.order <- match_arg(m.order, c("data", "random"))
    ord <- switch(m.order,
                  "random" = sample(1:n1, n1, replace = FALSE),
                  "data" = 1:n1)
  }

  for (r in 1:ratio) {
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

      #Get distances among eligible
      if (length(c.eligible) > 0) {
        #Compute PS differences if PS is distance or caliper
        if (!is.null(caliper) || !is.full.mahalanobis) {
          ps.diff <- abs(unname(d1[ord_i]) - d0[c.eligible])
        }

        if (is.null(mahcovs)) {
          #PS matching
          distances <- ps.diff
        }
        else {
          #MD matching
          distances <- sqrt(mahalanobis(mahcovs[c.eligible, ,drop = FALSE],
                                        mahcovs[lab1[ord_i],],
                                        cov = mahSigma_inv, inverted = TRUE))
        }
        # names(distances) <- c.eligible
      }

      if (!is.null(caliper) && length(c.eligible) > 0) {
        #PS caliper
        c.eligible <- c.eligible[ps.diff < cal]
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
              weights = weights.matrix(mm, treat, discarded))

  if (is.full.mahalanobis) class(res) <- c("matchit.mahalanobis", "matchit")
  else class(res) <- "matchit"

  return(res)
}

# MATCHIT method= subclass---------------------------------
#Needs updates
matchit2subclass <- function(treat, covs, data, distance, discarded,
                             ratio = 1, replace = FALSE, m.order = "largest",
                             caliper = NULL, mahvars = NULL, exact = NULL,
                             formula = NULL, verbose = FALSE,
                             is.full.mahalanobis,
                             ...) {

  if(verbose)
    cat("Subclassifying... \n")

  A <- list(...)
  subclass <- A[["subclass"]]
  sub.by <- A[["sub.by"]]

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

  if (is.null(sub.by)) sub.by <- "treat"
  else if (!is.vector(sub.by, "character") || length(sub.by) != 1){
    stop("sub.by must be one of \"treat\", \"control\", or \"all\".", call. = FALSE)
  }
  sub.by <- match_arg(sub.by, c("treat", "control", "all"))

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

  q <- switch(sub.by,
              "treat" = quantile(p1, probs = sprobs, na.rm=TRUE),
              "control" = quantile(p0, probs = sprobs, na.rm=TRUE),
              "all" = quantile(distance, probs = sprobs, na.rm=TRUE))

  ## Calculating Subclasses
  psclass <- setNames(rep(NA_integer_, n.obs), names(treat))
  psclass[!discarded] <- as.integer(findInterval(distance[!discarded], q, all.inside = TRUE))

  ## If any subclasses don't have members of a treatment group, fill them
  ## by "scooting" units from nearby subclasses until each subclass has a unit
  ## from each treatment group
  if (any(table(treat, psclass) == 0)) {
    psclass[!discarded] <- subclass_scoot(psclass[!discarded], treat[!discarded], distance[!discarded])
  }

  if (verbose) cat("Done\n")

  #warning for discrete data
  unique.classes <- unique(psclass[!is.na(psclass)], nmax = subclass)

  if (length(unique.classes) != subclass){
    warning("Due to discreteness in data, fewer subclasses generated than were requested.", call.=FALSE)
  }

  res <- list(subclass = psclass, q.cut = q,
              weights = weights.subclass(psclass, treat, estimand = "ATT"))

  class(res) <- c("matchit.subclass", "matchit")
  return(res)
}

