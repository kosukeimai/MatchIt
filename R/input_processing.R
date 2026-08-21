#Function to process inputs and throw warnings or errors if inputs are incompatible with methods
check.inputs <- function(mcall, method, distance, link, distance.options, exact, mahvars, antiexact,
                         caliper, discard, reestimate, s.weights, replace,
                         ratio, m.order, estimand, ...,
                         min.controls = NULL, max.controls = NULL) {

  null.method <- is_null(method)

  if (!null.method) {
    method <- arg::match_arg(method, c("exact", "cem", "nearest", "optimal", "full",
                                       "genetic", "subclass", "cardinality",
                                       "quick"))
  }

  ignored.inputs <- character(0L)
  error.inputs <- character(0L)

  .entered_arg <- function(mcall, i) {
    if (!hasName(mcall, i)) {
      return(FALSE)
    }

    i_ <- get0(i, envir = parent.frame(), inherits = FALSE)

    if (is_null(i_)) {
      return(FALSE)
    }

    !identical(i_, eval(formals(matchit)[[i]]))
  }

  if (null.method) {
    for (i in c("exact", "mahvars", "antiexact", "caliper", "std.caliper", "replace", "ratio", "min.controls", "max.controls", "m.order")) {
      if (.entered_arg(mcall, i)) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "exact") {
    for (i in c("distance", "link", "distance.options", "exact", "mahvars", "antiexact", "caliper", "std.caliper", "discard", "reestimate", "replace", "ratio", "min.controls", "max.controls", "m.order")) {
      if (.entered_arg(mcall, i)) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "cem") {
    for (i in c("distance", "link", "distance.options", "exact", "mahvars", "antiexact", "caliper", "std.caliper", "discard", "reestimate", "replace", "ratio", "min.controls", "max.controls")) {
      if (.entered_arg(mcall, i)) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "nearest") {
    if (rlang::is_string(distance) && distance %in% matchit_distances()) {
      for (e in c("mahvars", "reestimate")) {
        if (.entered_arg(mcall, e)) {
          error.inputs <- c(error.inputs, e)
        }
      }
    }
  }
  else if (method == "optimal") {
    if (rlang::is_string(distance) && distance %in% matchit_distances()) {
      for (e in c("mahvars", "reestimate")) {
        if (.entered_arg(mcall, e)) {
          error.inputs <- c(error.inputs, e)
        }
      }
    }

    for (i in c("replace", "caliper", "std.caliper", "m.order")) {
      if (.entered_arg(mcall, i)) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }

  }
  else if (method == "full") {
    if (rlang::is_string(distance) && distance %in% matchit_distances()) {
      for (e in c("mahvars", "reestimate")) {
        if (.entered_arg(mcall, e)) {
          error.inputs <- c(error.inputs, e)
        }
      }
    }

    for (i in c("replace", "ratio", "m.order")) {
      if (.entered_arg(mcall, i)) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "genetic") {
    if (rlang::is_string(distance) && distance %in% matchit_distances()) {
      for (e in c("mahvars", "reestimate")) {
        if (.entered_arg(mcall, e)) {
          error.inputs <- c(error.inputs, e)
        }
      }
    }
    for (i in c("min.controls", "max.controls")) {
      if (.entered_arg(mcall, i)) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "cardinality") {
    for (i in c("distance", "link", "distance.options", "antiexact", "caliper", "std.caliper", "reestimate", "replace", "min.controls", "m.order")) {
      if (.entered_arg(mcall, i)) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "subclass") {
    for (i in c("exact", "mahvars", "antiexact", "caliper", "std.caliper", "replace", "ratio", "min.controls", "max.controls", "m.order")) {
      if (.entered_arg(mcall, i)) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "quick") {
    if (rlang::is_string(distance) && distance %in% matchit_distances()) {
      for (e in c("mahvars", "reestimate")) {
        if (.entered_arg(mcall, e)) {
          error.inputs <- c(error.inputs, e)
        }
      }
    }

    for (i in c("replace", "ratio", "min.controls", "max.controls", "m.order", "antiexact")) {
      if (.entered_arg(mcall, i)) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }

  method_str <- if (null.method) "NULL" else dQuote(method, FALSE)

  if (is_not_null(ignored.inputs)) {
    arg::wrn("the {cli::qty(ignored.inputs)} argument{?s} {.arg {ignored.inputs}} {?is/are} not used with {.code method = {method_str}} and will be ignored")
  }

  if (is_not_null(error.inputs)) {
    distance_str <- dQuote(distance, FALSE)
    arg::err("the {cli::qty(error.inputs)} argument{?s} {.arg {error.inputs}} cannot be used with {.code method = {method_str}} and {.code distance = {distance_str}}")
  }

  ignored.inputs
}

#Check treatment for type, binary, missing, num. rows
check_treat <- function(treat = NULL, X = NULL) {

  if (is_null(treat)) {
    if (is_null(X) || is_null(attr(X, "treat"))) {
      return(NULL)
    }

    treat <- attr(X, "treat")
  }

  if (isTRUE(attr(treat, "checked"))) {
    return(treat)
  }

  if (!is.atomic(treat) || is_not_null(dim(treat))) {
    arg::err("the treatment must be a vector")
  }

  if (anyNA(treat)) {
    arg::err("missing values are not allowed in the treatment")
  }

  if (TRUE) {
    if (!has_n_unique(treat, 2L)) {
      arg::err("the treatment must be a binary variable")
    }

    treat <- binarize(treat) #make 0/1
  }
  else {
    if (has_n_unique(treat, 2L)) {
      treat <- {
        if (is.logical(treat) || all(as.character(treat) %in% c("0", "1")))
          factor(treat, levels = sort(unique(treat, nmax = 2)),
                 labels = c("control", "treated"), ordered = FALSE)
        else factor(treat, nmax = 2, ordered = FALSE)
      }

      # treat <- binarize(treat) #make 0/1

      attr(treat, "type") <- "binary"
      attr(treat, "treated") <- levels(treat)[2L]
      attr(treat, "ordered") <- FALSE
    }
    else {
      arg::err("the treatment must be a binary variable") #Remove to support multi

      if (!is.character(treat) && !is.factor(treat)) {
        arg::err("the treatment must be a factor variable if it takes on more than 2 unique values")
      }

      treat <- droplevels(as.factor(treat))

      attr(treat, "type") <- "multi"
      # attr(treat, "treated") <- levels(treat)[which.min(tabulateC(treat))]
      attr(treat, "ordered") <- is.ordered(treat)
    }
  }

  if (is_not_null(X) && length(treat) != nrow(X)) {
    arg::err("the treatment and covariates must have the same number of units")
  }

  attr(treat, "checked") <- TRUE

  treat
}

#Function to process distance and give warnings about new syntax
process.distance <- function(distance, method = NULL, treat) {
  if (is_null(distance)) {
    if (is_not_null(method) && !method %in% c("cem", "exact", "cardinality")) {
      arg::err("{.arg distance} cannot be {.val {list(NULL)}} with {.code method = {.str {method}}}")
    }

    return(distance)
  }

  if (rlang::is_string(distance)) {
    allowable.distances <- c(
      #Propensity score methods
      "glm", "cbps", "gam", "nnet", "rpart", "bart",
      "randomforest", "elasticnet", "lasso", "ridge", "gbm",
      #Distance matrices
      matchit_distances()
    )

    if (tolower(distance) %in% c("cauchit", "cloglog", "linear.cloglog", "linear.log", "linear.logit", "linear.probit",
                                 "linear.cauchit", "log", "probit")) {
      link <- tolower(distance)

      arg::wrn('{.code distance = {.str {distance}}} will be deprecated; please use {.code distance = "glm", link = {.str {link}}} in the future')

      distance <- "glm"
      attr(distance, "link") <- link
    }
    else if (tolower(distance) %in% tolower(c("GAMcloglog", "GAMlog", "GAMlogit", "GAMprobit"))) {
      link <- tolower(substr(distance, 4L, nchar(distance)))

      arg::wrn('{.code distance = {.str {distance}}} will be deprecated; please use {.code distance = "gam", link = {.str {link}}} in the future')

      distance <- "gam"
      attr(distance, "link") <- link
    }
    else if (tolower(distance) == "logit") {
      distance <- "glm"
      attr(distance, "link") <- "logit"
    }
    else if (tolower(distance) == "glmnet") {
      distance <- "elasticnet"
    }
    else if (!tolower(distance) %in% allowable.distances) {
      arg::err('the argument supplied to {.arg distance} is not an allowable value. See {.topic MatchIt::distance} for allowable options')
    }
    else if (is_not_null(method) && method == "subclass" && tolower(distance) %in% matchit_distances()) {
      arg::err('{.arg distance} cannot be {.val {distance}} with {.code method = "subclass"}')
    }
    else {
      distance <- tolower(distance)
    }

    return(distance)
  }

  if (!is.numeric(distance) || (is_not_null(dim(distance)) && length(dim(distance)) != 2)) {
    arg::err("{.arg distance} must be a string with the name of the distance measure to be used or a numeric vector or matrix containing distance measures")
  }

  if (is.matrix(distance) && (is_null(method) || !method %in% c("nearest", "optimal", "full"))) {
    method_str <- if (is_null(method)) "NULL" else dQuote(method, FALSE)
    arg::err("{.arg distance} cannot be supplied as a matrix with {.code method = {method_str}}")
  }

  if (is.matrix(distance)) {
    dim.distance <- dim(distance)

    if (all_equal_to(dim.distance, length(treat))) {
      if (is_not_null(rownames(distance))) {
        distance <- distance[names(treat), , drop = FALSE]
      }

      if (is_not_null(colnames(distance))) {
        distance <- distance[, names(treat), drop = FALSE]
      }

      distance <- distance[treat == 1, treat == 0, drop = FALSE]
    }
    else if (dim.distance[1L] == sum(treat == 1) &&
             dim.distance[2L] == sum(treat == 0)) {
      if (is_not_null(rownames(distance))) {
        distance <- distance[names(treat)[treat == 1], , drop = FALSE]
      }

      if (is_not_null(colnames(distance))) {
        distance <- distance[, names(treat)[treat == 0], drop = FALSE]
      }
    }
    else {
      arg::err("when supplied as a matrix, {.arg distance} must have dimensions NxN or N1xN0. See {.topic MatchIt::distance} for details")
    }
  }
  else if (length(distance) != length(treat)) {
    arg::err("{.arg distance} must be the same length as the dataset if specified as a numeric vector")
  }

  arg::arg_no_NA(distance)

  distance
}

#Function to check ratio is acceptable
process.ratio <- function(ratio, method = NULL, ..., min.controls = NULL, max.controls = NULL) {
  #Should be run after process.inputs() and ignored inputs set to NULL
  if (is_null(method)) {
    return(1)
  }

  ratio.null <- is_null(ratio)
  ratio.na <- !ratio.null && anyNA(ratio)

  if (method %in% c("nearest", "optimal")) {
    if (ratio.null) {
      ratio <- 1
    }

    arg::arg_number(ratio)
    arg::arg_gte(ratio, 1)

    if (is_null(max.controls)) {
      arg::arg_whole_number(ratio,
                            .msg = "{.arg ratio} must be a whole number when {.arg max.controls} is not specified")

      ratio <- round(ratio)
    }
    else {
      arg::arg_count(max.controls)

      if (ratio == 1) {
        arg::err("{.arg ratio} must be greater than 1 for variable ratio matching")
      }

      if (max.controls <= ratio) {
        arg::err("{.arg max.controls} must be greater than {.arg ratio} for variable ratio matching")
      }

      if (is_null(min.controls)) {
        min.controls <- 1
      }

      arg::arg_count(min.controls)

      if (min.controls < 1) {
        arg::err("{.arg min.controls} cannot be less than 1 for variable ratio matching")
      }

      if (min.controls >= ratio) {
        arg::err("{.arg min.controls} must be less than {.arg ratio} for variable ratio matching")
      }
    }
  }
  else if (method == "full") {
    if (is_null(max.controls)) {
      max.controls <- Inf
    }

    arg::arg_number(max.controls)
    arg::arg_gt(max.controls, 0)

    if (is_null(min.controls)) {
      min.controls <- 0
    }

    arg::arg_number(min.controls)
    arg::arg_gte(min.controls, 0)

    ratio <- 1 #Just to get min.controls and max.controls out
  }
  else if (method == "genetic") {
    if (ratio.null) {
      ratio <- 1
    }

    arg::arg_count(ratio)

    min.controls <- max.controls <- NULL
  }
  else if (method == "cardinality") {
    if (ratio.null) {
      ratio <- 1
    }
    else if (!ratio.na && (!is.numeric(ratio) || !identical(length(ratio), 1L) || ratio <= 0)) {
      arg::err("{.arg ratio} must be a single positive number or {.val {NA}}")
    }

    min.controls <- max.controls <- NULL
  }
  else {
    min.controls <- max.controls <- NULL
  }

  if (is_not_null(ratio)) {
    attr(ratio, "min.controls") <- min.controls
    attr(ratio, "max.controls") <- max.controls
  }

  ratio
}

#Function to check if caliper is okay and process it
process.caliper <- function(caliper = NULL, method = NULL, data = NULL, covs = NULL, mahcovs = NULL,
                            distance = NULL, discarded = NULL, std.caliper = TRUE) {

  #Check method; must be able to use a caliper
  #Check caliper names; if "" is one of them but distance = "mahal", throw error;
  #otherwise make sure variables exist in data or covs
  #Make sure no calipers are used on binary or factor variables (throw error if so)
  #Ignore calipers used on single-value variables or with caliper = NA or Inf
  #Export caliper.formula to add to covs
  #If std, export standardized versions

  #Check need for caliper
  if (is_null(caliper) || is_null(method) || !method %in% c("nearest", "genetic", "full", "quick")) {
    return(NULL)
  }

  #Check if form of caliper is okay
  arg::arg_numeric(caliper)

  #Check caliper names
  if (identical(length(caliper), 1L) && (is_null(names(caliper)) || identical(names(caliper), ""))) {
    names(caliper) <- ""
  }
  else if (is_null(names(caliper))) {
    arg::err("{.arg caliper} must be a named vector with names corresponding to the variables for which a caliper is to be applied")
  }
  else if (anyNA(names(caliper))) {
    arg::err("{.arg caliper} names cannot include {.val {NA}}")
  }
  else if (sum(!nzchar(names(caliper))) > 1L) {
    arg::err("no more than one entry in {.arg caliper} can have no name")
  }

  if (hasName(caliper, "") && is_null(distance)) {
    arg::err("all entries in {.arg caliper} must be named when {.arg distance} does not correspond to a propensity score")
  }

  #Check if caliper name is in available data
  cal.in.data <- setNames(names(caliper) %in% names(data), names(caliper))
  cal.in.covs <- setNames(names(caliper) %in% names(covs), names(caliper))
  cal.in.mahcovs <- setNames(names(caliper) %in% names(mahcovs), names(caliper))

  if (any(nzchar(names(caliper)) & !cal.in.covs & !cal.in.data)) {
    bad_vars <- names(caliper)[nzchar(names(caliper)) & !cal.in.data & !cal.in.covs & !cal.in.mahcovs]
    arg::err(c("All variables named in {.arg caliper} must be in {.arg data}.",
               "x" = "Variables not in {.arg data}: {.var {bad_vars}}"))
  }

  #Check std.caliper
  arg::arg_logical(std.caliper)

  if (length(std.caliper) == 1L) {
    std.caliper <- rep_with(std.caliper, caliper)
  }
  else if (length(std.caliper) == length(caliper)) {
    names(std.caliper) <- names(caliper)
  }
  else {
    arg::err("{.arg std.caliper} must have the same length as {.arg caliper}")
  }

  #Remove trivial calipers
  caliper <- caliper[is.finite(caliper)]

  if (is_null(caliper)) {
    return(NULL)
  }

  #Ensure no calipers on categorical variables
  cat.vars <- vapply(names(caliper), function(x) {
    v <- {
      if (!nzchar(x)) distance
      else if (cal.in.data[x]) data[[x]]
      else if (cal.in.covs[x]) covs[[x]]
      else mahcovs[[x]]
    }

    is.character(v) || is.factor(v)
  }, logical(1L))

  if (any(cat.vars)) {
    bad_vars <- ifelse(nzchar(names(caliper)), names(caliper), "<distance>")[cat.vars]
    arg::err(c("Calipers cannot be used with factor or character variables.",
               "x" = "Offending variables: {.var {bad_vars}}"))
  }

  #Process calipers according to std.caliper
  std.caliper <- std.caliper[names(std.caliper) %in% names(caliper)]

  arg::arg_no_NA(std.caliper)

  if (any(std.caliper)) {
    if (hasName(std.caliper, "") && isTRUE(std.caliper[!nzchar(names(std.caliper))]) && is.matrix(distance)) {
      arg::err("when {.arg distance} is supplied as a matrix and a caliper for it is specified, {.arg std.caliper} must be {.val {FALSE}} for the distance measure")
    }

    caliper[std.caliper] <- caliper[std.caliper] * vapply(names(caliper)[std.caliper], function(x) {
      if (!nzchar(x)) sd(distance[!discarded])
      else if (cal.in.data[x]) sd(data[[x]][!discarded])
      else if (cal.in.covs[x]) sd(covs[[x]][!discarded])
      else sd(mahcovs[[x]][!discarded])
    }, numeric(1L))
  }

  if (any(caliper < 0) && !method %in% c("nearest", "genetic", "full")) {
    arg::err("calipers cannot be negative with {.code method = {.str {method}}}")
  }

  #Add cal.formula
  if (any(nzchar(names(caliper)) & !cal.in.covs[names(caliper)] & !cal.in.mahcovs[names(caliper)])) {
    attr(caliper, "cal.formula") <- reformulate(names(caliper)[nzchar(names(caliper)) & !cal.in.covs[names(caliper)] & !cal.in.mahcovs[names(caliper)]])
  }

  caliper
}

#Function to process replace argument
process.replace <- function(replace, method = NULL, ..., reuse.max = NULL) {

  if (is_null(method)) {
    return(FALSE)
  }

  if (is_null(replace)) {
    replace <- FALSE
  }

  arg::arg_flag(replace)

  if (method %in% c("nearest")) {
    if (is_null(reuse.max)) {
      reuse.max <- if (replace) .Machine$integer.max else 1L
    }
    else {
      arg::arg_count(reuse.max)
      arg::arg_gte(reuse.max, 1)

      if (reuse.max > .Machine$integer.max) {
        reuse.max <- .Machine$integer.max
      }
    }

    replace <- reuse.max > 1L
    attr(replace, "reuse.max") <- as.integer(reuse.max)
  }

  replace
}

#Process variable input, e.g., to exact or mahvars, that accept a string or rhs formula
#Returns a model.frame object
process.variable.input <- function(x, data = NULL, n = rlang::caller_arg(x)) {

  if (is_null(x)) {
    return(NULL)
  }

  if (is.character(x)) {
    if (is_null(data) || !is.data.frame(data)) {
      arg::err("if {.arg {n}} is specified as strings, a data frame containing the named variables must be supplied to {.arg data}")
    }

    if (!all(hasName(data, x))) {
      arg::err(c("All names supplied to {.arg {n}} must be variables in {.arg data}.",
                 "x" = "Variables not in {.arg data}: {.var {setdiff(x, names(data))}}"))
    }

    x <- reformulate(x)
  }
  else if (rlang::is_formula(x)) {
    x <- update(terms(x, data = data), NULL ~ .)
  }
  else {
    arg::err("{.arg {n}} must be supplied as a character vector of names or a one-sided formula")
  }

  x_covs <- model.frame(x, data, na.action = "na.pass")

  if (anyNA(x_covs)) {
    arg::err("missing values are not allowed in the covariates named in {.arg {n}}")
  }

  x_covs
}
