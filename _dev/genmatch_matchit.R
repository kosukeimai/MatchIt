genmatch <- function(formula, data, treat, balance.covs.formula = NULL, reuse.max = 1, rat = 1, ord,
                     criterion = "smd.max", s.weights = NULL, cal.dist = NULL, dist = NULL, cal = NULL,
                     cal.covs = NULL, exact = NULL, anitexact = NULL, test = 4, ...) {
  #Argument checking
  if (missing(ord)) ord <- seq_len(sum(treat == 1))
  if (is.null(names(treat))) names(treat) <- rownames(data)
  if (length(rat) == 1) rat <- rep(rat, sum(treat == 1))
  if (reuse.max > sum(treat == 1)) reuse.max <- as.integer(sum(treat == 1))

  do_matching <- get(paste0("do_matching", test))

  #Initilize covariates for matching
  matching.covs <- transform_covariates(formula, data, treat = treat, method = "mahalanobis",
                                        s.weights = s.weights)

  #Initialize distance matrix restrictions
  matches_denied <- deny_match(treat, cal.dist, dist, cal, cal.covs, exact, anitexact)
  to_deny_matches <- length(matches_denied) > 0

  #Initilize covariates for balancing
  if (is.null(balance.covs.formula)) balance.covs.formula <- formula
  balance.covs <- get.covs.matrix(balance.covs.formula, data = data)

  #Initilize balance criterion
  bal.init <- initialize_balance(criterion = criterion, covs = balance.covs,
                                 treat = treat, s.weights = s.weights)

  #Function to do matching matching and produce balance statistic given covariate weights
  gen_fun <- function(w) {
    #Create distance matrix from weighted covariates
    distmat <- eucdist_internal(matching.covs %*% diag(w), treat)
    if (to_deny_matches) distmat[matches_denied] <- Inf

    #Perform matching
    matches <- do_matching(distmat, treat, ord, rat, reuse.max)

    #Extract matching weights
    match_weights <- weights.matrix(matches, treat)

    #Return balance statistic
    return(compute_balance(bal.init, match_weights))
  }

  #Set bounds on covariate weights; keep positive to limit search
  bounds <- matrix(c(0, 1000), nrow = ncol(matching.covs), ncol = 2,
                   byrow = TRUE)

  #Perform optimization
  opt <- rgenoud::genoud(gen_fun, nvars = ncol(matching.covs),
                         max = FALSE, Domains = bounds,
                         lexical = attr(bal.init, "lexical"),
                         gradient.check = FALSE, BFGS = FALSE,
                         hessian = FALSE, starting.values = rep(1, ncol(matching.covs)),
                         ...)

  #Once optimal weights found, create distance matrix
  distmat <- eucdist_internal(matching.covs %*% diag(opt$par), treat)
  distmat[matches_denied] <- Inf

  #Perform matching; output is a char match.matrix
  matches <- do_matching(distmat, treat, ord, rat, reuse.max)

  # browser()
  # microbenchmark::microbenchmark(
  #   `1` = do_matching1(distmat, treat, ord, rat, reuse.max),
  #   `1.5` = do_matching1.5(distmat, treat, ord, rat, reuse.max),
  #   `4` = do_matching4(distmat, treat, ord, rat, reuse.max),
  #   check = "equivalent", times = 100
  # )

  return(list(matches = matches, opt = opt, weights = weights.matrix(matches, treat)))
}

do_matching1 <- function(distmat, treat, ord = NULL, ratio, reuse.max) {
  #Matching w/ replacement quickly

  max_rat <- max(ratio)
  n1 <- nrow(distmat)
  mm <- matrix(NA_integer_, nrow = n1, ncol = max_rat,
               dimnames = list(rownames(distmat), NULL))
  num_matches <- rep(0, nrow(mm))

  for (r in 1:max_rat) {
    #get indices of minimum value
    m <- max.col(-distmat, "first")

    #inf: units where closest match is Inf away, i.e., not to be matched
    inf <- !is.finite(distmat[cbind(1:n1, m)])

    #Assign non-Inf matches
    mm[!inf, r] <- m[!inf]

    #Give matched pairs dist of Inf to prevent repeat matches
    distmat[cbind(which(!inf), m[!inf])] <- Inf

    #Increase number of matches for matched treated units by 1
    num_matches[!inf] <- num_matches[!inf] + 1

    #Treated units with enough matches get Inf distance for all control units
    if (r != max_rat) {
      distmat[num_matches == ratio,] <- Inf
    }
  }

  for (i in 1:ncol(mm)) {
    mm[!is.na(mm[,i]),i] <- which(treat == 0)[mm[!is.na(mm[,i]),i]]
  }

  mm <- nummm2charmm(mm, treat)

  return(mm)
}

do_matching1.5 <- function(distmat, treat, ord = NULL, ratio, reuse.max) {
  #Matching w/ replacement quickly
  #Implemented in C++ as do_matchingC2
  max_rat <- max(ratio)
  n1 <- nrow(distmat)
  mm <- matrix(NA_integer_, nrow = n1, ncol = max_rat,
               dimnames = list(rownames(distmat), NULL))
  num_matches_t <- rep(0, nrow(mm))
  num_matches_c <- rep(0, ncol(distmat))

  for (r in 1:max_rat) {
    for (i in which(num_matches_t < ratio)) {
      #get index of minimum value
      m <- which.min(distmat[i,])

      if (!is.finite(distmat[i, m])) next

      mm[i, r] <- m

      #Prevent re-matching
      distmat[i, m] <- Inf

      num_matches_t[i] <- num_matches_t[i] + 1
      num_matches_c[m] <- num_matches_c[m] + 1

      if (num_matches_c[m] == reuse.max) distmat[,m] <- Inf
    }
  }

  for (i in 1:ncol(mm)) {
    mm[!is.na(mm[,i]), i] <- which(treat == 0)[mm[!is.na(mm[,i]), i]]
  }

  mm <- nummm2charmm(mm, treat)

  return(mm)
}

do_matching4 <- function(distmat, treat, ord = NULL, ratio, reuse.max) {
  #Port of do_matching1.5
  mm <- do_matchingC(distmat, treat, ord, reuse.max, ratio)
  nummm2charmm(mm, treat)
}

#Returns vector of who is ineligible to be matched
deny_match <- function(treat, cal.dist = NULL, dist = NULL, cal = NULL, cal.covs = NULL, exact = NULL, antiexact = NULL) {
  d <- matrix(FALSE, nrow = sum(treat == 1), ncol = sum(treat == 0),
              dimnames = list(names(treat)[treat == 1], names(treat)[treat == 0]))

  if (!is.null(cal.dist)) {
    if (is.matrix(dist)) {
      d[d > cal.dist] <- TRUE
    }
    else {
      d[eucdist_internal(dist, treat) > cal.dist] <- TRUE
    }
  }

  if (!is.null(cal)) {
    for (i in seq_len(ncol(cal.covs))) {
      d[eucdist_internal(cal.covs[,i], treat) > cal[i]] <- TRUE
    }
  }

  if (!is.null(exact)) {
    d[outer(exact[treat==1], exact[treat==0], "!=")] <- TRUE
  }

  if (!is.null(antiexact)) {
    for (i in seq_len(ncol(cal.covs))) {
      d[outer(antiexact[treat==1, i], antiexact[treat==0, i], "==")] <- TRUE
    }
  }

  return(which(d))
}
