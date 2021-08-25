cardinality_matchit <- function(treat, X, estimand = "ATT", tols = .05, s.weights = NULL,
                                ratio = 1, solver = "glpk", time = 2*60, verbose = FALSE) {

  n <- length(treat)

  #Check inputs
  if (is.null(s.weights)) s.weights <- rep(1, n)
  else for (i in 0:1) s.weights[treat == i] <- s.weights[treat == i]/mean(s.weights[treat == i])

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
      0, 0       #slack coefs for each sample size (n1, n0)
    )

    #Constraint matrix
    target.means <- apply(X, 2, wm, w = s.weights)
    #One row per constraints, one column per coef
    C <- rbind(
      c(s.weights*(treat==1), -1, 0), #Num treated = n1
      c(s.weights*(treat!=1), 0, -1), #Num control = n0
      t(rbind((treat==1)*s.weights*X, -target.means-tols/2, 0)), #Balance constraints
      t(rbind((treat==1)*s.weights*X, -target.means+tols/2, 0)), #Means of each group must
      t(rbind((treat!=1)*s.weights*X, 0, -target.means-tols/2)), #be within tols/2 of full
      t(rbind((treat!=1)*s.weights*X, 0, -target.means+tols/2))  #sample means
    )

    #Constraint RHS
    #One per constraint
    Crhs <- c(0,
              0,
              rep(0, ncol(X)),
              rep(0, ncol(X)),
              rep(0, ncol(X)),
              rep(0, ncol(X))
    )

    #Constraint direction
    #One per constraint
    Cdir <- c("==",
              "==",
              rep("<", ncol(X)),
              rep(">", ncol(X)),
              rep("<", ncol(X)),
              rep(">", ncol(X))
    )

    #Coef types
    types <- c(rep("B", n), #Matching weights
               rep("C", 2)) #Slack coefs for matched group size (n1, n0)

    lower.bound <- c(rep(0, n),
                     rep(1, 2))

    #If ratio != 0, constrain n0 to be ratio*n1
    if (is.finite(ratio)) {
      C <- rbind(C, c(rep(0, n), ratio*1, -1))
      Crhs <- c(Crhs, 0)
      Cdir <- c(Cdir, "==")
    }
  }
  else if (match_type == "template_att") {
    #Find largest control group that matches treated group

    n0 <- sum(treat != 1)

    #Objective function: size of matched control group
    O <- c(
      rep(1, n0), #weights for each control unit
      0           #slack coef for size of control group
    )

    #Constraint matrix
    target.means <- apply(X[treat==1,,drop=FALSE], 2, wm, w = s.weights[treat==1])
    #One row per constraints, one column per coef
    C <- rbind(
      c(s.weights[treat==0], -1), #Num control = k
      t(rbind(s.weights[treat!=1]*X[treat!=1,,drop=FALSE], -target.means-tols)),
      t(rbind(s.weights[treat!=1]*X[treat!=1,,drop=FALSE], -target.means+tols))
    )

    #Constraint RHS
    #One per constraint
    Crhs <- c(
      0,
      rep(0, ncol(X)),
      rep(0, ncol(X))
    )

    #Constraint direction
    #One per constraint
    Cdir <- c(
      "==",
      rep("<", ncol(X)),
      rep(">", ncol(X))
    )

    #Coef types
    types <- c(rep("B", n0), #Matching weights
               rep("C", 1))  #Slack for num control matched

    lower.bound <- c(rep(0, n0),
                     rep(0, 1))
  }
  else if (match_type == "cardinality") {
    #True cardinality matching: find largest balanced sample

    #Objective function: total sample size
    O <- c(
      s.weights, #weight for each unit
      0          #coef for treated sample size (n1)
    )

    #Constraint matrix
    #One row per constraints, one column per coef
    C <- rbind(
      c(s.weights*(treat==1), -1),     #Num treated = n1
      c(s.weights*(treat!=1), -ratio), #Num control = ratio*n1
      t(rbind(((treat==1) - (treat!=1)/ratio)*s.weights*X, -tols)), #Balance constraints
      t(rbind(((treat==1) - (treat!=1)/ratio)*s.weights*X, tols))
    )

    #Constraint RHS
    #One per constraint
    Crhs <- c(
      0,
      0,
      rep(0, ncol(X)),
      rep(0, ncol(X))
    )

    #Constraint direction
    #One per constraint
    Cdir <- c("==",
              "==",
              rep("<", ncol(X)),
              rep(">", ncol(X))
    )

    #Coef types
    types <- c(rep("B", n), #Matching weights
               rep("C", 1)) #Slack coef for treated group size (n1)

    lower.bound <- c(rep(0, n),
                     rep(0, 1))

  }

  weights <- NULL

  opt.out <- dispatch_optimizer(solver = solver, obj = O, mat = C, dir = Cdir,
                                rhs = Crhs, types = types, max = TRUE, lb = lower.bound,
                                ub = NULL, time = time, verbose = verbose)

  cardinality_error_report(opt.out, solver)

  sol <- switch(solver, "glpk" = opt.out$solution,
                "symphony" = opt.out$solution,
                "gurobi" = opt.out$x)
  if (match_type %in% c("template_ate", "cardinality")) {
    weights <- sol[seq_len(n)]
  }
  else if (match_type %in% c("template_att")) {
    weights <- rep(1, n)
    weights[treat != 1] <- sol[seq_len(n0)]
  }

  #Make sure sum of weights in both groups is the same (important for exact matching)
  if (match_type == "template_att" && (is.na(ratio) || ratio != 1)) {
    weights[treat != 1] <- weights[treat != 1]*sum(weights[treat == 1])/sum(weights[treat != 1])
  }
  else {
    smaller.group <- c(0,1)[which.min(c(sum(weights[treat != 1]), sum(weights[treat == 1])))]
    weights[treat != smaller.group] <- weights[treat != smaller.group]*sum(weights[treat == smaller.group])/sum(weights[treat != smaller.group])
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