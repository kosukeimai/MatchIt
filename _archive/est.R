#Estimat effects and SEs
est_effect <- function(formula, m, type, measure, nboot = NULL, boot.type = "perc", robust.type = NULL, conf = .95, ...) {
  #formula: outcome model. Can include interactions because marginal effects will be computed after
  #m: matchit object, with call intact and with match.data available
  #type: type of SE/CI. Options: clusterrobust, robust, boot, blockboot
  #measure: effect measure. MD, RD, RR, OR, IRR
  #nboot: if bootstrap, how many
  #boot.type: if bootstrap, what kind of interval
  #robust.type: if robust, what type (HC1, etc.)
  #conf: confidence level

  A <- list(...)

  estimand <- m$estimand
  call <- m$call

  replace <- isTRUE(m$info$replace)
  method <- m$info$method

  #Extract matched data
  if (replace) {
    md <- do.call("get_matches", c(list(m), A[names(A) %in% names(formals(get_matches))]))
  }
  else {
    md <- do.call("match.data", c(list(m), A[names(A) %in% names(formals(match.data))]))
  }

  #Ensure treatment is 0/1
  treat.name <- as.character(m$formula[[2]])
  md[[treat.name]] <- binarize(md[[treat.name]])

  #Ensure treatment is in model formula
  tt.rhs <- delete.response(terms(formula))
  if (!treat.name %in% rownames(attr(tt.rhs, "factors"))) {
    stop(paste0("The treatment variable, ", treat.name, ", must be on the right-hand side of 'formula'."), call. = FALSE)
  }

  #Check if covariates are present
  covs.present <- any(rownames(attr(tt.rhs, "factors")) != treat.name)

  #Process outcome and get type
  outcome.mf <- model.frame(update(formula, . ~ 0), data = md)
  outcome <- model.response(outcome.mf)
  if (inherits(outcome, "Surv")) {
    out.type <- "survival"
    if (missing(measure)) measure <- "HR"
    else measure <- match_arg(measure, "HR")
  }
  else if ((is.factor(outcome) && nlevels(outcome) == 2) || length(unique(outcome)) == 2) {
    out.type <- "binary"
    if (missing(measure)) measure <- "RD"
    else measure <- match_arg(measure, c("RD", "RR", "OR"))
  }
  else {
    out.type <- "continuous"
    if (missing(measure)) measure <- "MD"
    else measure <- match_arg(measure, "MD")
  }

  #Process SE type
  if (missing(type)) {
    if (method %in% c("subclass", "exact", "cem")) type <- "robust"
    else type <- "clusterrobust"
  }
  else {
    type <- match_arg(type, c("clusterrobust", "robust", "boot", "blockboot"))
  }

  if (type %in% c("boot", "blockboot")) {
    if (missing(nboot) || is.null(nboot)) stop(paste0("An argument to 'nboot' must be supplied when type = \"", type, "\"."), call. = FALSE)
  }
  else nboot <- NULL

  #Check compatibility and assign glm functions
  if (out.type == "continuous") {
    fit.fun <- function(f, d, w) {
      environment(f) <- environment()
      lm(f, data = d, weights = w)
    }
  }
  else if (out.type == "binary") {
    if (is.null(nboot) && covs.present) stop("Bootstrapping must be used when including covariates in the model for a binary outcome.\n\tSet 'type' to \"blockboot\" or \"boot\" and supply an argument to 'nboot'.", call. = FALSE)
    fit.fun <- function(f, d, w) {
      environment(f) <- environment()
      if (is.null(nboot)) {
        if (measure == "RD") lm(f, data = d, weights = w) #Cheung (2007)
        else if (measure == "RR") glm(f, data = d, weights = w, family = quasipoisson(link = "log")) #Chen et al (2018)
        else if (measure == "OR") glm(f, data = d, weights = w, family = quasibinomial(link = "logit"))
      }
      else {
        glm(f, data = d, weights = w, family = quasibinomial(link = "logit"))
      }
    }
  }
  else if (out.type == "survival") {
    if (covs.present) stop("Covariates cannot be included in the outcome model with survival outcomes.", call. = FALSE)
    fit.fun <- function(f, d, w) {
      environment(f) <- environment()
      if (replace && is.null(nboot)) {
        #Austin & Cafri's (2020) SE estimator
        fs <- survival::coxph(f, data = d, robust = TRUE, weights = w,
                              cluster = data[[attr(d, "subclass")]])
        Vs <- fs$var
        ks <- nlevels(data[[attr(d, "subclass")]])

        fi <- survival::coxph(f, data = d, robust = TRUE, weights = w,
                              cluster = data[[attr(d, "id")]])
        Vi <- fi$var
        ki <- length(unique(data[[attr(d, "id")]]))

        fc <- survival::coxph(f, data = d, robust = TRUE, weights = w)
        Vc <- fc$var
        kc <- nrow(d)

        #Compute the variance
        V <- (ks/(ks-1))*Vs + (ki/(ki-1))*Vi - (kc/(kc-1))*Vc

        #Sneak it back into the fit object
        fc$var <- V
      }
      else {
        fc <- survival::coxph(f, data = d, robust = is.null(nboot), weights = w,
                              cluster = if (type == "clusterrobust") d[[attr(d, "subclass")]])
      }
      return(fc)
    }
  }

  if (is.null(nboot)) {
    w1 <- md[[attr(md, "weights")]]
    fit <- fit.fun(formula, md, w1)

    if (out.type == "continuous") {

      beta <- coef(fit)

      if (type == "robust") sigma <- sandwich::vcovHC(fit, type = robust.type)
      else if (type == "clusterrobust") {
        if (replace) sigma <- sandwich::vcovCL(fit, cluster = md[c(attr(md, "subclass"), attr(md, "id"))], type = robust.type)
        else sigma <- sandwich::vcovCL(fit, cluster = md[[attr(md, "subclass")]], type = robust.type)
      }

      if (estimand == "ATT") {
        md <- md[md[[treat.name]] == 1,]
      } else if (estimand == "ATC") {
        md <- md[md[[treat.name]] == 0,]
      }

      md[[treat.name]] <- 1
      alpha1 <- md[[attr(md, "weights")]] %*% model.matrix(formula, md) / sum(md[[attr(md, "weights")]])
      md[[treat.name]] <- 0
      alpha0 <- md[[attr(md, "weights")]] %*% model.matrix(formula, md) / sum(md[[attr(md, "weights")]])

      alpha_tau <- alpha1 - alpha0

      alpha <- rbind(alpha1, alpha0, alpha_tau)

      est <- alpha %*% beta
      vc <- alpha %*% sigma %*% t(alpha)
      se <- sqrt(diag(vc))

      df <- fit$df.residual

    }
    else if (out.type == "binary") {
      beta <- coef(fit)

      if (replace) sigma <- sandwich::vcovCL(fit, cluster = md[c(attr(md, "subclass"), attr(md, "id"))], type = robust.type)
      else sigma <- sandwich::vcovCL(fit, cluster = md[[attr(md, "subclass")]], type = robust.type)

      alpha <- matrix(c(1,1,
                        1,0,
                        0,1), nrow = 3, byrow = TRUE)

      est <- alpha %*% beta
      vc <- alpha %*% sigma %*% t(alpha)
      se <- sqrt(diag(vc))

      df <- Inf

    }
    else if (out.type == "survival") {
      if (!is.null(robust.type)) warning("'robust.type' is ignored with survival outcomes.", call. = FALSE, immediate. = TRUE)
      est <- coef(fit)[treat.name]
      vc <- vcov(fit)[treat.name, treat.name, drop = FALSE]
      se <- sqrt(diag(vc))

      df <- Inf
    }

    #Compute effect measures on link scale
    l.effects <- cbind(est,
                     est + qt(.5*(1-conf), df) * se,
                     est + qt(1 - .5*(1-conf), df) * se,
                     se,
                     tval <- est/se,
                     2 * pt(abs(tval), df, lower.tail = FALSE))
  }
  else {
    if (!is.null(robust.type)) warning("'robust.type' is ignored when using bootstrapping.", call. = FALSE, immediate. = TRUE)

    if (missing(boot.type) || is.null(boot.type)) {
      warning("'boot.type' not specificed. Using \"perc\".", call. = FALSE, immediate. = TRUE)
      boot.type <- "perc"
    }
    else {
      if (identical(boot.type, "stud")) stop("Studentized bootstrap CIs are not allowed. Set 'boot.type' to something other than \"stud\".", call. = FALSE)
      boot.type <- match_arg(boot.type, c("norm", "basic", "perc", "bca"))
    }

    fun1 <- function(x) switch(measure, MD = x, RD = x, RR = log(x), OR = log(x/(1-x)))

    if (type == "boot") {

      if (!is.null(m$model)) {
        env <- attributes(terms(m$model))$.Environment
      } else {
        env <- parent.frame()
      }
      data <- eval(call$data, envir = env)

      if (boot.type == "bca" && nboot < nrow(data)) {
        stop(paste0("'nboot' cannot be less than the number of rows of the original dataset (", nrow(data),
                    ") when type = \"boot\" and boot.type = \"bca\". Increase 'nboot' or change 'boot.type'."), call. = FALSE)
      }

      estfun <- function(data, i) {
        boot_call <- call
        boot_call[["data"]] <- data[i,]

        m_boot <- eval(boot_call)

        md_boot <- do.call("match.data", c(list(m_boot), A[names(A) %in% names(formals(match.data))]))

        fit <- fit.fun(formula, md_boot, md_boot[[attr(md_boot, "weights")]])

        if (out.type == "survival") {
          return(coef(fit)[treat.name])
        }
        else {
          if (estimand == "ATT") {
            md_boot <- md_boot[md_boot[[treat.name]] == 1,]
          } else if (estimand == "ATC") {
            md_boot <- md_boot[md_boot[[treat.name]] == 0,]
          }

          md_boot[[treat.name]] <- 1
          M1 <- weighted.mean(predict(fit, newdata = md_boot, type = "response"), md_boot[[attr(md_boot, "weights")]])

          md_boot[[treat.name]] <- 0
          M0 <- weighted.mean(predict(fit, newdata = md_boot, type = "response"), md_boot[[attr(md_boot, "weights")]])

          return(c(fun1(M1), fun1(M0), fun1(M1) - fun1(M0)))
        }
      }

      boot.out <- do.call(boot::boot, c(list(data, estfun, R = nboot), A[names(A) %in% names(formals(boot::boot))]))
    }
    else if (type == "blockboot") {

      pairs <- levels(md[[attr(md, "subclass")]])
      pairs.num <- seq_along(pairs)

      if (boot.type == "bca" && nboot < length(pairs)) {
        stop(paste0("'nboot' cannot be less than the number of matched pairs/subclass (", length(pairs),
                    ") when type = \"blockboot\" and boot.type = \"bca\". Increase 'nboot' or change 'boot.type'."), call. = FALSE)
      }

      estfun <- function(data, i) {

        #Compute number of times each pair is present
        numreps <- setNames(tabulate(data[i], length(pairs)), pairs)

        #For each pair p, copy corresponding md row indices numreps[p] times
        ids <- unlist(lapply(pairs[numreps > 0], function(p) rep(which(md[[attr(md, "subclass")]] == p), numreps[p])))

        #Subset md with block bootstrapped ids
        md_boot <- md[ids,]

        fit <- fit.fun(formula, md_boot, md_boot[[attr(md_boot, "weights")]])

        if (out.type == "survival") {
          return(coef(fit)[treat.name])
        }
        else {
          if (estimand == "ATT") {
            md_boot <- md_boot[md_boot[[treat.name]] == 1,]
          } else if (estimand == "ATC") {
            md_boot <- md_boot[md_boot[[treat.name]] == 0,]
          }

          md_boot[[treat.name]] <- 1
          M1 <- weighted.mean(predict(fit, newdata = md_boot, type = "response"), md_boot[[attr(md_boot, "weights")]])

          md_boot[[treat.name]] <- 0
          M0 <- weighted.mean(predict(fit, newdata = md_boot, type = "response"), md_boot[[attr(md_boot, "weights")]])

          return(c(fun1(M1), fun1(M0), fun1(M1) - fun1(M0)))
        }
      }

      boot.out <- do.call(boot::boot, c(list(pairs.num, estfun, R = nboot), A[names(A) %in% names(formals(boot::boot))]))
    }

    est <- boot.out$t0
    vc <- cov(boot.out$t)
    se <- sqrt(diag(vc))

    ci <- t(vapply(seq_along(est), function(x) {
      bc <- boot::boot.ci(boot.out, conf = conf, type = boot.type, index = x)
      drop(bc[[4]][,-seq_len(length(bc[[4]]) - 2)])
    }, numeric(2)))

    l.effects <- cbind(est, ci, se, tval <- est/se, 2 * pnorm(abs(tval), lower.tail = FALSE))

    df <- Inf
  }

  dimnames(l.effects) <- list(switch(measure,
                                   MD = c("E[Y1]", "E[Y0]", "E[Y1] - E[Y0]"),
                                   RD = c("P(Y1=1)", "P(Y0=1)", "RD"),
                                   RR = c("log(P(Y1=1))", "log(P(Y0=1))", "log(RR)"),
                                   OR = c("logit(P(Y1=1))", "logit(P(Y0=1))", "log(OR)"),
                                   HR = c("log(HR)")),
                            c("Estimate", "CI lower", "CI upper", "Std. Error",
                              if (is.finite(df)) "t value" else "z value",
                              if (is.finite(df)) "Pr(>|t|)" else "Pr(>|z|)"))

  effects <- switch(measure,
                    MD = l.effects[,1:3],
                    RD = l.effects[,1:3],
                    RR = exp(l.effects[,1:3]),
                    OR = rbind(exp(l.effects[1:2,1:3])/(1 + exp(l.effects[1:2,1:3])), exp(l.effects[3,1:3])),
                    HR = exp(l.effects[,1:3, drop = FALSE]))

  if (!is.null(effects)) {
    dimnames(effects) <- list(switch(measure,
                                     MD = c("E[Y1]", "E[Y0]", "E[Y1] - E[Y0]"),
                                     HR = "HR",
                                     c("P(Y1=1)", "P(Y0=1)", measure)),
                              c("Estimate",
                                "CI lower", "CI upper"))
  }

  dimnames(vc) <- list(rownames(l.effects), rownames(l.effects))

  out <- list(effects = effects, l.effects = l.effects, vcov = vc)
  attr(out, "info") <- list(out.type = out.type, covs.present = covs.present,
                            replace = replace, type = type,
                            boot.type = boot.type, robust.type = robust.type,
                            estimand = estimand, df = df)
  class(out) <- "eff_matchit"

  out
}

print.eff_matchit <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat(paste("\nMarginal Effect Estimates:\n\n", sep = ""))
  printCoefmat(x$effects, digits = digits, has.Pvalue = FALSE, right = TRUE, cs.ind = 1:3,
               tst.ind = NULL, ...)

  cat(paste("\nMarginal Effect Estimates (Link Scale):\n\n", sep = ""))
  printCoefmat(x$l.effects, digits = digits, cs.ind = 1:4, tst.ind = which(endsWith(colnames(x$effects), " value")),
               has.Pvalue = TRUE, right = TRUE, signif.stars = FALSE, ...)
  cat("\n")

  #Print info
  # info <- attr(x, "info")
  # if (info) 0
  # "SEs and CIs estimated on the link scale using HC1 cluster-robust variance"
  # "SEs estimated on link scale using block bootstrap; percentile CIs used"
  invisible(x)
}

coef.eff_matchit <- function(object, link = TRUE, ...) {
  if (link) setNames(object$l.effects[,"Estimate"], rownames(object$l.effects))
  else setNames(object$effects[,"Estimate"], rownames(object$effects))
}

vcov.eff_matchit <- function(object, ...) {
  object$vcov
}