cardinality_matchit <- function(treat, X, estimand = "ATT", tols = .05, s.weights = NULL,
                                ratio = 1, focal = NULL, tvals = NULL, solver = "glpk", time = 2*60, verbose = FALSE) {

  n <- length(treat)
  if (is.null(tvals)) tvals <- unique(treat)
  nt <- length(tvals)

  #Check inputs
  if (is.null(s.weights)) s.weights <- rep(1, n)
  else for (i in tvals) s.weights[treat == i] <- s.weights[treat == i]/mean(s.weights[treat == i])

  if (is.null(focal)) focal <- max(tvals)

  if (length(time) != 1 || !is.numeric(time) || time <= 0) stop("'time' must be a positive number.", call. = FALSE)

  solver <- match_arg(solver, c("glpk", "symphony", "gurobi"))
  check.package(switch(solver, glpk = "Rglpk", symphony = "Rsymphony", gurobi = "gurobi"))

  #Select match type
  if (estimand == "ATE") match_type <- "template_ate"
  else if (!is.finite(ratio)) match_type <- "template_att"
  else match_type <- "cardinality"

  #Set objective and constraints
  if (match_type == "template_ate") {
    #Find largest sample that matches full sample

    #Objective function: total sample size
    O <- c(
      s.weights, #weight for each unit
      rep(0, nt)       #slack coefs for each sample size (n1, n0)
    )

    #Constraint matrix
    target.means <- apply(X, 2, wm, w = s.weights)

    C <- matrix(0, nrow = nt * (1 + 2*ncol(X)), ncol = length(O))
    Crhs <- rep(0, nrow(C))
    Cdir <- rep("==", nrow(C))

    for (i in seq_len(nt)) {
      #Num in group i = ni
      C[i, seq_len(n)] <- s.weights * (treat == tvals[i])
      C[i, n + i] <- -1

      #Cov means must be less than target.means+tols/2
      r1 <- nt + (i - 1)*2*ncol(X) + 1:ncol(X)
      C[r1, seq_len(n)] <- t((treat==tvals[i])*s.weights*X)
      C[r1, n + i] <- -target.means-tols/2
      Cdir[r1] <- "<"

      #Cov means must be greater than target.means-tols/2
      r2 <- r1 + ncol(X)
      C[r2, seq_len(n)] <- t((treat==tvals[i])*s.weights*X)
      C[r2, n + i] <- -target.means+tols/2
      Cdir[r2] <- ">"
    }

    #If ratio != 0, constrain n0 to be ratio*n1
    if (nt == 2L && is.finite(ratio)) {
      C_ratio <- c(rep(0, n), rep(-1, nt))
      C_ratio[n + which(tvals == focal)] <- ratio
      C <- rbind(C, C_ratio)
      Crhs <- c(Crhs, 0)
      Cdir <- c(Cdir, "==")
    }

    #Coef types
    types <- c(rep("B", n), #Matching weights
               rep("C", nt)) #Slack coefs for matched group size

    lower.bound <- c(rep(0, n),
                     rep(1, nt))
    upper.bound <- NULL
  }
  else if (match_type == "template_att") {
    #Find largest control group that matches treated group

    nonf <- which(treat != focal)
    n0 <- length(nonf)
    tvals_ <- setdiff(tvals, focal)

    #Objective function: size of matched control group
    O <- c(
      rep(1, n0), #weights for each non-focal unit
      rep(0, nt - 1)  #slack coef for size of non-focal groups
    )

    #Constraint matrix
    target.means <- apply(X[treat==focal,,drop=FALSE], 2, wm, w = s.weights[treat==focal])
    #One row per constraint, one column per coef

    C <- matrix(0, nrow = (nt - 1) * (1 + 2*ncol(X)), ncol = length(O))
    Crhs <- rep(0, nrow(C))
    Cdir <- rep("==", nrow(C))

    for (i in seq_len(nt - 1)) {
      #Num in group i = ni
      C[i, seq_len(n0)] <- s.weights[nonf] * (treat[nonf] == tvals_[i])
      C[i, n0 + i] <- -1

      #Cov means must be less than target.means+tols
      r1 <- nt - 1 + (i - 1)*2*ncol(X) + 1:ncol(X)
      C[r1, seq_len(n0)] <- t((treat[nonf]==tvals_[i])*s.weights[nonf]*X[nonf,,drop = FALSE])
      C[r1, n0 + i] <- -target.means-tols
      Cdir[r1] <- "<"

      #Cov means must be greater than target.means-tols
      r2 <- r1 + ncol(X)
      C[r2, seq_len(n0)] <- t((treat[nonf]==tvals_[i])*s.weights[nonf]*X[nonf,,drop = FALSE])
      C[r2, n0 + i] <- -target.means+tols
      Cdir[r2] <- ">"
    }

    #Coef types
    types <- c(rep("B", n0), #Matching weights
               rep("C", nt - 1))  #Slack for num control matched

    lower.bound <- c(rep(0, n0),
                     rep(0, nt - 1))
    upper.bound <- NULL
  }
  else if (match_type == "cardinality") {
    #True cardinality matching: find largest balanced sample
    if (nt > 2) ratio <- 1

    #Objective function: total sample size
    O <- c(
      s.weights, #weight for each unit
      0          #coef for treated sample size (n1)
    )

    #Constraint matrix
    t_combs <- combn(tvals, 2, simplify = FALSE)

    C <- matrix(0, nrow = nt + 2*ncol(X)*length(t_combs), ncol = length(O))
    Crhs <- rep(0, nrow(C))
    Cdir <- rep("==", nrow(C))

    for (i in seq_len(nt)) {
      #Num in group i = ni
      C[i, seq_len(n)] <- s.weights * (treat == tvals[i])
      C[i, n + 1] <- if (tvals[i] == focal) -1 else -ratio
    }

    for (j in seq_along(t_combs)) {
      t_comb <- t_combs[[j]]
      if (t_comb[2] == focal) t_comb <- rev(t_comb)

      r1 <- nt + (j - 1)*2*ncol(X) + 1:ncol(X)
      C[r1, seq_len(n)] <- t(((treat==t_comb[1]) - (treat==t_comb[2])/ratio)*s.weights*X)
      C[r1, n + 1] <- -tols
      Cdir[r1] <- "<"

      r2 <- r1 + ncol(X)
      C[r2, seq_len(n)] <- t(((treat==t_comb[1]) - (treat==t_comb[2])/ratio)*s.weights*X)
      C[r2, n + 1] <- tols
      Cdir[r2] <- ">"
    }

    #Coef types
    types <- c(rep("B", n), #Matching weights
               rep("C", 1)) #Slack coef for treated group size (n1)

    lower.bound <- c(rep(0, n),
                     rep(0, 1))
    upper.bound <- NULL
  }

  weights <- NULL

  opt.out <- dispatch_optimizer(solver = solver, obj = O, mat = C, dir = Cdir,
                                rhs = Crhs, types = types, max = TRUE, lb = lower.bound,
                                ub = upper.bound, time = time, verbose = verbose)

  cardinality_error_report(opt.out, solver)

  sol <- switch(solver, "glpk" = opt.out$solution,
                "symphony" = opt.out$solution,
                "gurobi" = opt.out$x)
  if (match_type %in% c("template_ate", "cardinality")) {
    weights <- round(sol[seq_len(n)])
  }
  else if (match_type %in% c("template_att")) {
    weights <- rep(1, n)
    weights[treat != focal] <- round(sol[seq_len(n0)])
  }

  #Make sure sum of weights in both groups is the same (important for exact matching)
  if (match_type == "template_att" && (is.na(ratio) || ratio != 1)) {
    for (t in setdiff(tvals, focal)) {
      weights[treat == t] <- weights[treat == t]*sum(weights[treat == focal])/sum(weights[treat == t])
    }
  }
  else {
    smallest.group <- tvals[which.min(vapply(tvals, function(t) sum(treat == t), numeric(1L)))]
    for (t in setdiff(tvals, smallest.group)) {
      weights[treat == t] <- weights[treat == t]*sum(weights[treat == smallest.group])/sum(weights[treat == t])
    }
  }

  return(list(weights = weights, opt.out = opt.out))
}

cardinality_error_report <- function(out, solver) {
  if (solver == "glpk") {
    if (out$status == 1) {
      if (all(out$solution == 0)) {
        stop("The optimization problem may be infeasible. Try increasing the value of 'tols'.\nSee ?method_cardinality for additional details.", call. = FALSE)
      }
      else {
        warning("The optimizer failed to find an optimal solution in the time alotted. The returned solution may not be optimal.\nSee ?method_cardinality for additional details.", call. = FALSE)
      }
    }
  }
  else if (solver == "symphony") {
    if (names(out$status) %in% c("TM_TIME_LIMIT_EXCEEDED") && !all(out$solution == 0) && all(out$solution <= 1)) {
      warning("The optimizer failed to find an optimal solution in the time alotted. The returned solution may not be optimal.", call. = FALSE)
    }
    else if (names(out$status) != "TM_OPTIMAL_SOLUTION_FOUND") {
      stop("The optimizer failed to find an optimal solution in the time alotted. The optimization problem may be infeasible. Try increasing the value of 'tols'.\nSee ?method_cardinality for additional details.", call. = FALSE)
    }
  }
  else if (solver == "gurobi") {
    if (out$status %in% c("TIME_LIMIT", "SUBOPTIMAL") && !all(out$x == 0)) {
      warning("The optimizer failed to find an optimal solution in the time alotted. The returned solution may not be optimal.\nSee ?method_cardinality for additional details.", call. = FALSE)
    }
    else if (out$status %in% c("INFEASIBLE", "INF_OR_UNBD", "NUMERIC") || all(out$x == 0)) {
      stop("The optimization problem may be infeasible. Try increasing the value of 'tols'.\nSee ?method_cardinality for additional details.", call. = FALSE)
    }
  }
}

dispatch_optimizer <- function(solver = "glpk", obj, mat, dir, rhs, types, max = TRUE, lb = NULL, ub = NULL, time = NULL, verbose = FALSE) {
  if (solver == "glpk") {
    dir[dir == "="] <- "=="
    opt.out <- Rglpk::Rglpk_solve_LP(obj = obj, mat = mat, dir = dir, rhs = rhs, max = max,
                                     types = types,
                                     # bounds = list(lower = lb, upper = ub), #Spurious warning when using bounds
                                     control = list(tm_limit = time*1000, verbose = verbose))
  }
  else if (solver == "symphony") {
    dir[dir == "<"] <- "<="
    dir[dir == ">"] <- ">="
    dir[dir == "="] <- "=="
    opt.out <- Rsymphony::Rsymphony_solve_LP(obj = obj, mat = mat, dir = dir, rhs = rhs, max = TRUE,
                                             types = types, verbosity = verbose - 2,
                                             # bounds = list(lower = lb, upper = ub), #Spurious warning when using bounds
                                             time_limit = time)
  }
  else if (solver == "gurobi") {
    dir[dir == "<="] <- "<"
    dir[dir == ">="] <- ">"
    dir[dir == "=="] <- "="
    opt.out <- gurobi::gurobi(list(A = mat, obj = obj, sense = dir, rhs = rhs, vtype = types,
                                   modelsense = "max", lb = lb, ub = ub),
                              params = list(OutputFlag = as.integer(verbose), TimeLimit = time))
  }

  opt.out
}
