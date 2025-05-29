#Function to process inputs and throw warnings or errors if inputs are incompatible with methods
check.inputs <- function(mcall, method, distance, link, distance.options, exact, mahvars, antiexact,
                         caliper, discard, reestimate, s.weights, replace,
                         ratio, m.order, estimand, ...,
                         min.controls = NULL, max.controls = NULL) {

  null.method <- is_null(method)

  if (null.method) {
    method <- "NULL"
  }
  else {
    method <- match_arg(method, c("exact", "cem", "nearest", "optimal", "full",
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
    if (is.character(distance) && distance %in% matchit_distances()) {
      for (e in c("mahvars", "reestimate")) {
        if (.entered_arg(mcall, e)) {
          error.inputs <- c(error.inputs, e)
        }
      }
    }
  }
  else if (method == "optimal") {
    if (is.character(distance) && distance %in% matchit_distances()) {
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
    if (is.character(distance) && distance %in% matchit_distances()) {
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
    if (is.character(distance) && distance %in% matchit_distances()) {
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
    if (is.character(distance) && distance %in% matchit_distances()) {
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

  if (is_not_null(ignored.inputs)) {
    .wrn(sprintf("the argument%%s %s %%r not used with `method = %s` and will be ignored",
                 word_list(ignored.inputs, quotes = "`"),
                 add_quotes(method, quotes = !null.method)),
         n = length(ignored.inputs))
  }

  if (is_not_null(error.inputs)) {
    .err(sprintf("the argument%%s %s %%r not used with `method = %s` and `distance = %s`",
                 word_list(error.inputs, quotes = "`"),
                 add_quotes(method, quotes = !null.method),
                 add_quotes(distance)),
         n = length(error.inputs))
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
    .err("the treatment must be a vector")
  }

  if (anyNA(treat)) {
    .err("missing values are not allowed in the treatment")
  }

  if (TRUE) {
    if (!has_n_unique(treat, 2L)) {
      .err("the treatment must be a binary variable")
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
      .err("the treatment must be a binary variable") #Remove to support multi

      if (!chk::vld_character_or_factor(treat)) {
        .err("the treatment must be a factor variable if it takes on more than 2 unique values")
      }

      treat <- droplevels(as.factor(treat))

      attr(treat, "type") <- "multi"
      # attr(treat, "treated") <- levels(treat)[which.min(tabulateC(treat))]
      attr(treat, "ordered") <- is.ordered(treat)
    }
  }

  if (is_not_null(X) && length(treat) != nrow(X)) {
    .err("the treatment and covariates must have the same number of units")
  }

  attr(treat, "checked") <- TRUE

  treat
}

#Function to process distance and give warnings about new syntax
process.distance <- function(distance, method = NULL, treat) {
  if (is_null(distance)) {
    if (is_not_null(method) && !method %in% c("cem", "exact", "cardinality")) {
      .err(sprintf("`distance` cannot be `NULL` with `method = %s`",
                   add_quotes(method)))
    }

    return(distance)
  }

  if (chk::vld_string(distance)) {
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

      .wrn(sprintf('`distance = "%s"` will be deprecated; please use `distance = "glm", link = "%s"` in the future',
                   distance, link))

      distance <- "glm"
      attr(distance, "link") <- link
    }
    else if (tolower(distance) %in% tolower(c("GAMcloglog", "GAMlog", "GAMlogit", "GAMprobit"))) {
      link <- tolower(substr(distance, 4L, nchar(distance)))

      .wrn(sprintf('`distance = "%s"` will be deprecated; please use `distance = "gam", link = "%s"` in the future',
                   distance, link))

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
      .err('the argument supplied to `distance` is not an allowable value. See `help("distance", package = "MatchIt")` for allowable options')
    }
    else if (is_not_null(method) && method == "subclass" && tolower(distance) %in% matchit_distances()) {
      .err(sprintf('`distance` cannot be %s with `method = "subclass"`',
                   add_quotes(distance)))
    }
    else {
      distance <- tolower(distance)
    }

    return(distance)
  }

  if (!is.numeric(distance) || (is_not_null(dim(distance)) && length(dim(distance)) != 2)) {
    .err("`distance` must be a string with the name of the distance measure to be used or a numeric vector or matrix containing distance measures")
  }

  if (is.matrix(distance) && (is_null(method) || !method %in% c("nearest", "optimal", "full"))) {
    .err(sprintf("`distance` cannot be supplied as a matrix with `method = %s`",
                 add_quotes(method, quotes = is_not_null(method))))
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
      .err("when supplied as a matrix, `distance` must have dimensions NxN or N1xN0. See `help(\"distance\")` for details")
    }
  }
  else if (length(distance) != length(treat)) {
    .err("`distance` must be the same length as the dataset if specified as a numeric vector")
  }

  chk::chk_not_any_na(distance)

  distance
}

#Function to check ratio is acceptable
process.ratio <- function(ratio, method = NULL, ..., min.controls = NULL, max.controls = NULL) {
  #Should be run after process.inputs() and ignored inputs set to NULL
  ratio.null <- is_null(ratio)
  ratio.na <- !ratio.null && anyNA(ratio)

  if (is_null(method)) {
    return(1)
  }

  if (method %in% c("nearest", "optimal")) {
    if (ratio.null) {
      ratio <- 1
    }
    else {
      chk::chk_number(ratio)
      chk::chk_gte(ratio, 1)
    }

    if (is_null(max.controls)) {
      if (!chk::vld_whole_number(ratio)) {
        .err("`ratio` must be a whole number when `max.controls` is not specified")
      }

      ratio <- round(ratio)
    }
    else {
      chk::chk_count(max.controls)

      if (ratio == 1) {
        .err("`ratio` must be greater than 1 for variable ratio matching")
      }

      if (max.controls <= ratio) {
        .err("`max.controls` must be greater than `ratio` for variable ratio matching")
      }

      if (is_null(min.controls)) {
        min.controls <- 1
      }
      else {
        chk::chk_count(min.controls)
      }

      if (min.controls < 1) {
        .err("`min.controls` cannot be less than 1 for variable ratio matching")
      }

      if (min.controls >= ratio) {
        .err("`min.controls` must be less than `ratio` for variable ratio matching")
      }
    }
  }
  else if (method == "full") {
    if (is_null(max.controls)) {
      max.controls <- Inf
    }
    else {
      chk::chk_number(max.controls)
      chk::chk_gt(max.controls, 0)
    }

    if (is_null(min.controls)) {
      min.controls <- 0
    }
    else {
      chk::chk_number(min.controls)
      chk::chk_gt(min.controls, 0)
    }

    ratio <- 1 #Just to get min.controls and max.controls out
  }
  else if (method == "genetic") {
    if (ratio.null) {
      ratio <- 1
    }
    else {
      chk::chk_count(ratio)
    }

    min.controls <- max.controls <- NULL
  }
  else if (method == "cardinality") {
    if (ratio.null) {
      ratio <- 1
    }
    else if (!ratio.na && (!chk::vld_number(ratio) || !chk::vld_gte(ratio, 0))) {
      .err("`ratio` must be a single positive number or `NA`")
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
  if (!is.atomic(caliper) || !is.numeric(caliper)) {
    .err("`caliper` must be a numeric vector")
  }

  #Check caliper names
  if (length(caliper) == 1L && (is_null(names(caliper)) || identical(names(caliper), ""))) {
    names(caliper) <- ""
  }
  else if (is_null(names(caliper))) {
    .err("`caliper` must be a named vector with names corresponding to the variables for which a caliper is to be applied")
  }
  else if (anyNA(names(caliper))) {
    .err("`caliper` names cannot include `NA`")
  }
  else if (sum(!nzchar(names(caliper))) > 1L) {
    .err("no more than one entry in `caliper` can have no name")
  }

  if (hasName(caliper, "") && is_null(distance)) {
    .err("all entries in `caliper` must be named when `distance` does not correspond to a propensity score")
  }

  #Check if caliper name is in available data
  cal.in.data <- setNames(names(caliper) %in% names(data), names(caliper))
  cal.in.covs <- setNames(names(caliper) %in% names(covs), names(caliper))
  cal.in.mahcovs <- setNames(names(caliper) %in% names(mahcovs), names(caliper))
  if (any(nzchar(names(caliper)) & !cal.in.covs & !cal.in.data)) {
    .err(paste0("All variables named in `caliper` must be in `data`. Variables not in `data`:\n\t",
                toString(names(caliper)[nzchar(names(caliper)) & !cal.in.data & !cal.in.covs & !cal.in.mahcovs])),
         tidy = FALSE)
  }

  #Check std.caliper
  chk::chk_logical(std.caliper)
  if (length(std.caliper) == 1L) {
    std.caliper <- rep_with(std.caliper, caliper)
  }
  else if (length(std.caliper) == length(caliper)) {
    names(std.caliper) <- names(caliper)
  }
  else {
    .err("`std.caliper` must be the same length as `caliper`")
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

    chk::vld_character_or_factor(v)
  }, logical(1L))

  if (any(cat.vars)) {
    .err(paste0("Calipers cannot be used with factor or character variables. Offending variables:\n\t",
                toString(ifelse(nzchar(names(caliper)), names(caliper), "<distance>")[cat.vars])),
         tidy = FALSE)
  }

  #Process calipers according to std.caliper
  std.caliper <- std.caliper[names(std.caliper) %in% names(caliper)]
  chk::chk_not_any_na(std.caliper)

  if (any(std.caliper)) {
    if (hasName(std.caliper, "") && isTRUE(std.caliper[!nzchar(names(std.caliper))]) && is.matrix(distance)) {
      .err("when `distance` is supplied as a matrix and a caliper for it is specified, `std.caliper` must be `FALSE` for the distance measure")
    }

    caliper[std.caliper] <- caliper[std.caliper] * vapply(names(caliper)[std.caliper], function(x) {
      if (!nzchar(x)) sd(distance[!discarded])
      else if (cal.in.data[x]) sd(data[[x]][!discarded])
      else if (cal.in.covs[x]) sd(covs[[x]][!discarded])
      else sd(mahcovs[[x]][!discarded])
    }, numeric(1L))
  }

  if (any(caliper < 0) && !method %in% c("nearest", "genetic", "full")) {
    .err(sprintf("calipers cannot be negative with `method = %s`",
                 add_quotes(method)))
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

  chk::chk_flag(replace)

  if (method %in% c("nearest")) {
    if (is_null(reuse.max)) {
      reuse.max <- if (replace) .Machine$integer.max else 1L
    }
    else {
      chk::chk_count(reuse.max)
      chk::chk_gte(reuse.max, 1)

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
process.variable.input <- function(x, data = NULL) {
  n <- deparse1(substitute(x))

  if (is_null(x)) {
    return(NULL)
  }

  if (is.character(x)) {
    if (is_null(data) || !is.data.frame(data)) {
      .err(sprintf("if `%s` is specified as strings, a data frame containing the named variables must be supplied to `data`",
                   n))
    }

    if (!all(hasName(data, x))) {
      .err(sprintf("All names supplied to `%s` must be variables in `data`. Variables not in `data`:\n\t%s", n,
                   toString(add_quotes(setdiff(x, names(data))))), tidy = FALSE)
    }

    x <- reformulate(x)
  }
  else if (rlang::is_formula(x)) {
    x <- update(terms(x, data = data), NULL ~ .)
  }
  else {
    .err(sprintf("`%s` must be supplied as a character vector of names or a one-sided formula.", n))
  }

  x_covs <- model.frame(x, data, na.action = "na.pass")

  if (anyNA(x_covs)) {
    .err(sprintf("missing values are not allowed in the covariates named in `%s`", n))
  }

  x_covs
}
