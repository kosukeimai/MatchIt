#Function to process inputs and throw warnings or errors if inputs are incompatible with methods
check.inputs <- function(mcall, method, distance, exact, mahvars, antiexact,
                         caliper, discard, reestimate, s.weights, replace,
                         ratio, m.order, estimand, ...,
                         min.controls = NULL, max.controls = NULL) {

  null.method <- is.null(method)
  if (null.method) {
    method <- "NULL"
  }
  else {
    method <- match_arg(method, c("exact", "cem", "nearest", "optimal", "full", "genetic", "subclass", "cardinality",
                                  "quick"))
  }

  ignored.inputs <- character(0)
  error.inputs <- character(0)
  if (null.method) {
    for (i in c("exact", "mahvars", "antiexact", "caliper", "std.caliper", "replace", "ratio", "min.controls", "max.controls", "m.order")) {
      if (i %in% names(mcall) && !is.null(i_ <- get0(i, inherits = FALSE)) && !identical(i_, formals(matchit)[[i]])) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "exact") {
    for (i in c("distance", "exact", "mahvars", "antiexact", "caliper", "std.caliper", "discard", "reestimate", "replace", "ratio", "min.controls", "max.controls", "m.order")) {
      if (i %in% names(mcall) && !is.null(i_ <- get0(i, inherits = FALSE)) && !identical(i_, formals(matchit)[[i]])) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "cem") {
    for (i in c("distance", "exact", "mahvars", "antiexact", "caliper", "std.caliper", "discard", "reestimate", "replace", "ratio", "min.controls", "max.controls", "m.order")) {
      if (i %in% names(mcall) && !is.null(i_ <- get0(i, inherits = FALSE)) && !identical(i_, formals(matchit)[[i]])) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "nearest") {
    if (is.character(distance) && distance %in% matchit_distances()) {
      for (e in c("mahvars", "reestimate")) {
        if (e %in% names(mcall) && !is.null(e_ <- get0(e, inherits = FALSE)) && !identical(e_, formals(matchit)[[e]])) {
          error.inputs <- c(error.inputs, e)
        }
      }
    }
  }
  else if (method == "optimal") {
    if (is.character(distance) && distance %in% matchit_distances()) {
      for (e in c("mahvars", "reestimate")) {
        if (e %in% names(mcall) && !is.null(e_ <- get0(e, inherits = FALSE)) && !identical(e_, formals(matchit)[[e]])) {
          error.inputs <- c(error.inputs, e)
        }
      }
    }

    for (i in c("replace", "caliper", "std.caliper", "m.order")) {
      if (i %in% names(mcall) && !is.null(i_ <- get0(i, inherits = FALSE)) && !identical(i_, formals(matchit)[[i]])) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }

  }
  else if (method == "full") {
    if (is.character(distance) && distance %in% matchit_distances()) {
      for (e in c("mahvars", "reestimate")) {
        if (e %in% names(mcall) && !is.null(e_ <- get0(e, inherits = FALSE)) && !identical(e_, formals(matchit)[[e]])) {
          error.inputs <- c(error.inputs, e)
        }
      }
    }

    for (i in c("replace", "ratio", "m.order")) {
      if (i %in% names(mcall) && !is.null(i_ <- get0(i, inherits = FALSE)) && !identical(i_, formals(matchit)[[i]])) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "genetic") {
    if (is.character(distance) && distance %in% matchit_distances()) {
      for (e in c("mahvars", "reestimate")) {
        if (e %in% names(mcall) && !is.null(e_ <- get0(e, inherits = FALSE)) && !identical(e_, formals(matchit)[[e]])) {
          error.inputs <- c(error.inputs, e)
        }
      }
    }
    for (i in c("min.controls", "max.controls")) {
      if (i %in% names(mcall) && !is.null(i_ <- get0(i, inherits = FALSE)) && !identical(i_, formals(matchit)[[i]])) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "cardinality") {
    for (i in c("distance", "antiexact", "caliper", "std.caliper", "reestimate", "replace", "min.controls", "m.order")) {
      if (i %in% names(mcall) && !is.null(i_ <- get0(i, inherits = FALSE)) && !identical(i_, formals(matchit)[[i]])) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "subclass") {
    for (i in c("exact", "mahvars", "antiexact", "caliper", "std.caliper", "replace", "ratio", "min.controls", "max.controls", "m.order")) {
      if (i %in% names(mcall) && !is.null(i_ <- get0(i, inherits = FALSE)) && !identical(i_, formals(matchit)[[i]])) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "quick") {
    if (is.character(distance) && distance %in% matchit_distances()) {
      for (e in c("mahvars", "reestimate")) {
        if (e %in% names(mcall) && !is.null(e_ <- get0(e, inherits = FALSE)) && !identical(e_, formals(matchit)[[e]])) {
          error.inputs <- c(error.inputs, e)
        }
      }
    }

    for (i in c("replace", "ratio", "min.controls", "max.controls", "m.order", "antiexact")) {
      if (i %in% names(mcall) && !is.null(i_ <- get0(i, inherits = FALSE)) && !identical(i_, formals(matchit)[[i]])) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }

  if (length(ignored.inputs) > 0) .wrn(sprintf("the %s %s not used with `method = %s` and will be ignored",
                                                  ngettext(length(ignored.inputs), "argument", "arguments"),
                                                  word_list(ignored.inputs, quotes = 1, is.are = TRUE),
                                                  add_quotes(method, quotes = !null.method)))
  if (length(error.inputs) > 0) .err(sprintf("the %s %s not used with `method = %s` and `distance = \"%s\"`",
                                             ngettext(length(error.inputs), "argument", "arguments"),
                                             word_list(error.inputs, quotes = 1, is.are = TRUE),
                                             add_quotes(method, quotes = !null.method),
                                             distance))
  ignored.inputs
}

#Check treatment for type, binary, missing, num. rows
check_treat <- function(treat = NULL, X = NULL) {

  if (is.null(treat)) {
    if (is.null(X) || is.null(attr(X, "treat"))) return(NULL)
    treat <- attr(X, "treat")
  }
  if (isTRUE(attr(treat, "checked"))) return(treat)

  if (!is.atomic(treat) || !is.null(dim(treat))) {
    .err("the treatment must be a vector")
  }

  if (anyNA(treat)) .err("missing values are not allowed in the treatment")
  if (length(unique(treat)) != 2) .err("the treatment must be a binary variable")
  if (!is.null(X) && length(treat) != nrow(X)) .err("the treatment and covariates must have the same number of units")

  treat <- binarize(treat) #make 0/1
  attr(treat, "checked") <- TRUE
  treat
}

#Function to process distance and give warnings about new syntax
process.distance <- function(distance, method = NULL, treat) {
  if (is.null(distance)) {
    if (!is.null(method) && !method %in% c("cem", "exact", "cardinality")) {
      .err(sprintf("`distance` cannot be `NULL` with `method = \"%s\"`",
                   method))
    }
    return(distance)
  }

  if (is.character(distance) && length(distance) == 1) {
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
      .wrn(sprintf("`distance = \"%s\"` will be deprecated; please use `distance = \"glm\", link = \"%s\"` in the future",
                      distance, link))
      distance <- "glm"
      attr(distance, "link") <- link
    }
    else if (tolower(distance) %in% tolower(c("GAMcloglog", "GAMlog", "GAMlogit", "GAMprobit"))) {
      link <- tolower(substr(distance, 4, nchar(distance)))
      .wrn(sprintf("`distance = \"%s\"` will be deprecated; please use `distance = \"gam\", link = \"%s\"` in the future",
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
      .err("the argument supplied to `distance` is not an allowable value. See `help(\"distance\")` for allowable options")
    }
    else if (!is.null(method) && method == "subclass" && tolower(distance) %in% matchit_distances()) {
      .err(sprintf("`distance` cannot be %s with `method = \"subclass\"`",
                   add_quotes(distance)))
    }
    else {
      distance <- tolower(distance)
    }

  }
  else if (!is.numeric(distance) || (!is.null(dim(distance)) && length(dim(distance)) != 2)) {
    .err("`distance` must be a string with the name of the distance measure to be used or a numeric vector or matrix containing distance measures")
  }
  else if (is.matrix(distance) && (is.null(method) || !method %in% c("nearest", "optimal", "full"))) {
    .err(sprintf("`distance` cannot be supplied as a matrix with `method = %s`",
                 add_quotes(method, quotes = !is.null(method))))
  }

  if (is.numeric(distance)) {
    if (is.matrix(distance)) {
      dim.distance <- dim(distance)
      if (all(dim.distance == length(treat))) {
        if (!is.null(rownames(distance))) distance <- distance[names(treat),, drop = FALSE]
        if (!is.null(colnames(distance))) distance <- distance[,names(treat), drop = FALSE]
        distance <- distance[treat == 1, treat == 0, drop = FALSE]
      }
      else if (all(dim.distance == c(sum(treat==1), sum(treat==0)))) {
        if (!is.null(rownames(distance))) distance <- distance[names(treat)[treat == 1],, drop = FALSE]
        if (!is.null(colnames(distance))) distance <- distance[,names(treat)[treat == 0], drop = FALSE]
      }
      else {
        .err("when supplied as a matrix, `distance` must have dimensions NxN or N1xN0. See `help(\"distance\")` for details")
      }
    }
    else {
      if (length(distance) != length(treat)) {
        .err("`distance` must be the same length as the dataset if specified as a numeric vector")
      }
    }

    chk::chk_not_any_na(distance)
  }

  distance
}

#Function to check ratio is acceptable
process.ratio <- function(ratio, method = NULL, ..., min.controls = NULL, max.controls = NULL) {
  #Should be run after process.inputs() and ignored inputs set to NULL
  ratio.null <- length(ratio) == 0
  ratio.na <- !ratio.null && anyNA(ratio)

  if (is.null(method)) return(1)
  if (method %in% c("nearest", "optimal")) {
    if (ratio.null) ratio <- 1
    else if (ratio.na) .err("`ratio` cannot be `NA`")
    else if (!is.atomic(ratio) || !is.numeric(ratio) || length(ratio) > 1 || ratio < 1) {
      .err("`ratio` must be a single number greater than or equal to 1")
    }

    if (is.null(max.controls)) {
      if (!chk::vld_whole_number(ratio)) {
        .err("`ratio` must be a whole number when `max.controls` is not specified")
      }
      ratio <- round(ratio)
    }
    else if (anyNA(max.controls) || !is.atomic(max.controls) || !is.numeric(max.controls) || length(max.controls) > 1) {
      .err("`max.controls` must be a single positive number")
    }
    else {
      if (ratio <= 1) .err("`ratio` must be greater than 1 for variable ratio matching")

      max.controls <- ceiling(max.controls)
      if (max.controls <= ratio) .err("`max.controls` must be greater than `ratio` for variable ratio matching")

      if (is.null(min.controls)) min.controls <- 1
      else if (anyNA(max.controls) || !is.atomic(max.controls) || !is.numeric(max.controls) || length(max.controls) > 1) {
        .err("`max.controls` must be a single positive number")
      }
      else min.controls <- floor(min.controls)

      if (min.controls < 1) .err("`min.controls` cannot be less than 1 for variable ratio matching")
      if (min.controls >= ratio) .err("`min.controls` must be less than `ratio` for variable ratio matching")
    }
  }
  else if (method == "full") {
    if (is.null(max.controls)) max.controls <- Inf
    else if ((anyNA(max.controls) || !is.atomic(max.controls) || !is.numeric(max.controls) || length(max.controls) > 1)) {
      .err("`max.controls` must be a single positive number")
    }

    if (is.null(min.controls)) min.controls <- 0
    else if ((anyNA(min.controls) || !is.atomic(min.controls) || !is.numeric(min.controls) || length(min.controls) > 1)) {
      .err("`min.controls` must be a single positive number")
    }

    ratio <- 1 #Just to get min.controls and max.controls out
  }
  else if (method == "genetic") {
    if (ratio.null) ratio <- 1
    else if (ratio.na) .err("`ratio` cannot be `NA`")
    else if (!is.atomic(ratio) || !is.numeric(ratio) || length(ratio) > 1 || ratio < 1 ||
             !chk::vld_whole_number(ratio)) {
      .err("`ratio` must be a single whole number greater than or equal to 1")
    }
    ratio <- round(ratio)

    min.controls <- max.controls <- NULL
  }
  else if (method == "cardinality") {
    if (ratio.null) ratio <- 1
    else if (!ratio.na && (!is.atomic(ratio) || !is.numeric(ratio) || length(ratio) > 1 || ratio < 0)) {
      .err("`ratio` must be a single positive number or `NA`")
    }

    min.controls <- max.controls <- NULL
  }
  else {
    min.controls <- max.controls <- NULL
  }

  if (!is.null(ratio)) {
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
  if (length(caliper) == 0 || is.null(method) || !method %in% c("nearest", "genetic", "full", "quick")) return(NULL)

  #Check if form of caliper is okay
  if (!is.atomic(caliper) || !is.numeric(caliper)) .err("`caliper` must be a numeric vector")

  #Check caliper names
  if (length(caliper) == 1 && (is.null(names(caliper)) || identical(names(caliper), ""))) {
    names(caliper) <- ""
  }
  else if (is.null(names(caliper))) {
    .err("`caliper` must be a named vector with names corresponding to the variables for which a caliper is to be applied")
  }
  else if (anyNA(names(caliper))) {
    .err("`caliper` names cannot include `NA`")
  }
  else if (sum(names(caliper) == "") > 1) {
    .err("no more than one entry in `caliper` can have no name")
  }

  if (any(names(caliper) == "") && is.null(distance)) {
    .err("all entries in `caliper` must be named when `distance` does not correspond to a propensity score")
  }

  #Check if caliper name is in available data
  cal.in.data <- setNames(names(caliper) %in% names(data), names(caliper))
  cal.in.covs <- setNames(names(caliper) %in% names(covs), names(caliper))
  cal.in.mahcovs <- setNames(names(caliper) %in% names(mahcovs), names(caliper))
  if (any(names(caliper) != "" & !cal.in.covs & !cal.in.data)) {
    .err(paste0("All variables named in `caliper` must be in `data`. Variables not in `data`:\n\t",
                paste0(names(caliper)[names(caliper) != "" & !cal.in.data & !cal.in.covs & !cal.in.mahcovs], collapse = ", ")), tidy = FALSE)
  }

  #Check std.caliper
  chk::chk_logical(std.caliper)
  if (length(std.caliper) == 1) {
    std.caliper <- setNames(rep.int(std.caliper, length(caliper)), names(caliper))
  }
  else if (length(std.caliper) != length(caliper)) {
    .err("`std.caliper` must be the same length as `caliper`")
  }
  else names(std.caliper) <- names(caliper)

  #Remove trivial calipers
  caliper <- caliper[is.finite(caliper)]

  num.unique <- vapply(names(caliper), function(x) {
    if (x == "") var <- distance
    else if (cal.in.data[x]) var <- data[[x]]
    else if (cal.in.covs[x]) var <- covs[[x]]
    else var <- mahcovs[[x]]

    length(unique(var))
  }, integer(1L))

  caliper <- caliper[num.unique > 1]

  if (length(caliper) == 0) return(NULL)

  #Ensure no calipers on categorical variables
  cat.vars <- vapply(names(caliper), function(x) {
    if (num.unique[names(num.unique) == x] == 2) return(TRUE)

    if (x == "") var <- distance
    else if (cal.in.data[x]) var <- data[[x]]
    else if (cal.in.covs[x]) var <- covs[[x]]
    else var <- mahcovs[[x]]

    is.factor(var) || is.character(var)
  }, logical(1L))

  if (any(cat.vars)) {
    .err(paste0("Calipers cannot be used with binary, factor, or character variables. Offending variables:\n\t",
                paste0(ifelse(names(caliper) == "", "<distance>", names(caliper))[cat.vars], collapse = ", ")),
         tidy = FALSE)
  }

  #Process calipers according to std.caliper
  std.caliper <- std.caliper[names(std.caliper) %in% names(caliper)]
  chk::chk_not_any_na(std.caliper)

  if (any(std.caliper)) {
    if ("" %in% names(std.caliper) && isTRUE(std.caliper[names(std.caliper) == ""]) && is.matrix(distance)) {
      .err("when `distance` is supplied as a matrix and a caliper for it is specified, `std.caliper` must be `FALSE` for the distance measure")
    }
    caliper[std.caliper] <- caliper[std.caliper] * vapply(names(caliper)[std.caliper], function(x) {
      if (x == "") sd(distance[!discarded])
      else if (cal.in.data[x]) sd(data[[x]][!discarded])
      else if (cal.in.covs[x]) sd(covs[[x]][!discarded])
      else sd(mahcovs[[x]][!discarded])
    }, numeric(1L))
  }

  #Add cal.formula
  if (any(names(caliper) != "" & !cal.in.covs[names(caliper)] & !cal.in.mahcovs[names(caliper)])) {
    attr(caliper, "cal.formula") <- reformulate(names(caliper)[names(caliper) != "" & !cal.in.covs[names(caliper)] & !cal.in.mahcovs[names(caliper)]])
  }

  abs(caliper)
}

#Function to process replace argument
process.replace <- function(replace, method = NULL, ..., reuse.max = NULL) {

  if (is.null(method)) return(FALSE)

  if (is.null(replace)) replace <- FALSE
  chk::chk_flag(replace)

  if (method %in% c("nearest")) {
    if (is.null(reuse.max)) {
      if (replace) reuse.max <- .Machine$integer.max
      else reuse.max <- 1L
    }
    else if (length(reuse.max) == 1 && is.numeric(reuse.max) &&
             (!is.finite(reuse.max) || reuse.max > .Machine$integer.max) &&
             !anyNA(reuse.max)) {
      reuse.max <- .Machine$integer.max
    }
    else if (abs(reuse.max - round(reuse.max)) > 1e-8 || length(reuse.max) != 1 ||
             anyNA(reuse.max) || reuse.max < 1) {
      .err("`reuse.max` must be a positive integer of length 1")
    }

    replace <- reuse.max != 1L
    attr(replace, "reuse.max") <- as.integer(reuse.max)
  }

  replace
}

#Process variable input, e.g., to exact or mahvars, that accept a string or rhs formula
#Returns a model.frame object
process.variable.input <- function(x, data = NULL) {
  n <- deparse1(substitute(x))

  if (is.null(x)) return(NULL)

  if (is.character(x)) {
    if (is.null(data) || !is.data.frame(data)) {
      .err(sprintf("if `%s` is specified as strings, a data frame containing the named variables must be supplied to `data`",
                   n))
    }
    if (!all(x %in% names(data))) {
      .err(sprintf("All names supplied to `%s` must be variables in `data`. Variables not in `data`:\n\t%s", n,
                   paste(add_quotes(setdiff(x, names(data))), collapse = ", ")), tidy = FALSE)
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