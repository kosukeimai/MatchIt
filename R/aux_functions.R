#Auxiliary functions; some from WeightIt

#Function to process inputs and throw warnings or errors if inputs are incompatible with methods
check.inputs <- function(mcall, method, distance, exact, mahvars, antiexact, caliper, discard, reestimate, s.weights, replace, ratio, min.controls = NULL, max.controls = NULL, m.order, estimand, ...) {

  null.method <- is.null(method)
  if (null.method) {
    method <- "NULL"
  }
  else {
    method <- match_arg(method, c("exact", "cem", "nearest", "optimal", "full", "genetic", "subclass", "cardinality"))
  }

  ignored.inputs <- character(0)
  error.inputs <- character(0)
  if (null.method) {
    for (i in c("exact", "mahvars", "antiexact", "caliper", "std.caliper", "replace", "ratio", "min.controls", "max.controls", "m.order")) {
      if (i %in% names(mcall) && !is.null(get0(i, inherits = FALSE))) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "exact") {
    for (i in c("distance", "exact", "mahvars", "antiexact", "caliper", "std.caliper", "discard", "reestimate", "replace", "ratio", "min.controls", "max.controls", "m.order")) {
      if (i %in% names(mcall) && !is.null(get0(i, inherits = FALSE))) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "cem") {
    for (i in c("distance", "exact", "mahvars", "antiexact", "caliper", "std.caliper", "discard", "reestimate", "replace", "ratio", "min.controls", "max.controls", "m.order")) {
      if (i %in% names(mcall) && !is.null(get0(i, inherits = FALSE))) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "nearest") {
    if (is.character(distance) && distance == "mahalanobis") {
      for (e in c("mahvars", "reestimate")) {
        if (e %in% names(mcall) && !is.null(get0(e, inherits = FALSE))) {
          error.inputs <- c(error.inputs, e)
        }
      }
    }
  }
  else if (method == "optimal") {
    if (is.character(distance) && distance == "mahalanobis") {
      for (e in c("mahvars", "reestimate")) {
        if (e %in% names(mcall) && !is.null(get0(e, inherits = FALSE))) {
          error.inputs <- c(error.inputs, e)
        }
      }
    }

    for (i in c("replace", "caliper", "std.caliper", "m.order")) {
      if (i %in% names(mcall) && !is.null(get0(i, inherits = FALSE))) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }

  }
  else if (method == "full") {
    if (is.character(distance) && distance == "mahalanobis") {
      for (e in c("mahvars", "reestimate")) {
        if (e %in% names(mcall) && !is.null(get0(e, inherits = FALSE))) {
          error.inputs <- c(error.inputs, e)
        }
      }
    }

    for (i in c("replace", "ratio", "m.order")) {
      if (i %in% names(mcall) && !is.null(get0(i, inherits = FALSE))) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "genetic") {
    if (is.character(distance) && distance == "mahalanobis") {
      for (e in c("mahvars", "reestimate")) {
        if (e %in% names(mcall) && !is.null(get0(e, inherits = FALSE))) {
          error.inputs <- c(error.inputs, e)
        }
      }
    }
    for (i in c("min.controls", "max.controls")) {
      if (i %in% names(mcall) && !is.null(get0(i, inherits = FALSE))) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "cardinality") {
    for (i in c("distance", "mahvars", "antiexact", "caliper", "std.caliper", "reestimate", "replace", "min.controls", "m.order")) {
      if (i %in% names(mcall) && !is.null(get0(i, inherits = FALSE))) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "subclass") {
    if (is.character(distance) && distance == "mahalanobis") {
      stop("distance = \"mahalanobis\" is not compatible with subclassification.", call. = FALSE)
    }

    for (i in c("exact", "mahvars", "antiexact", "caliper", "std.caliper", "replace", "ratio", "min.controls", "max.controls", "m.order")) {
      if (i %in% names(mcall) && !is.null(get0(i, inherits = FALSE))) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }

  if (length(ignored.inputs) > 0) warning(paste0(ngettext(length(ignored.inputs), "The argument ", "The arguments "),
                                                 word_list(ignored.inputs, quotes = 1, is.are = TRUE),
                                                 " not used with method = ", add_quotes(method, quotes = !null.method),
                                                 " and will be ignored."),
                                          call. = FALSE, immediate. = TRUE)
  if (length(error.inputs) > 0) stop(paste0(ngettext(length(error.inputs), "The argument ", "The arguments "),
                                            word_list(error.inputs, quotes = 1, is.are = TRUE),
                                            " not allowed with method = ", add_quotes(method, quotes = !null.method),
                                            " and distance = \"", distance, "\"."),
                                     call. = FALSE)
  return(ignored.inputs)
}

#Function to process distance and give warnings about new syntax
process.distance <- function(distance, method, treat) {
  if (is.null(distance) && !is.null(method)) stop(paste0("'distance' cannot be NULL with method = \"", method, "\"."), call. = FALSE)
  else if (is.character(distance) && length(distance) == 1) {
    allowable.distances <- c("glm", "cbps", "gam", "mahalanobis", "nnet", "rpart", "bart", "randomforest", "glmnet", "gbm")

    if (tolower(distance) %in% c("cauchit", "cloglog", "linear.cloglog", "linear.log", "linear.logit", "linear.probit",
                        "linear.cauchit", "log", "probit")) {
      warning(paste0("'distance = \"", distance, "\"' will be deprecated; please use 'distance = \"glm\", link = \"", distance, "\"' in the future."), call. = FALSE, immediate. = TRUE)
      link <- distance
      distance <- "glm"
      attr(distance, "link") <- link
    }
    else if (tolower(distance) %in% tolower(c("GAMcloglog", "GAMlog", "GAMlogit", "GAMprobit"))) {
      link <- sub("GAM", "", distance)
      warning(paste0("'distance = \"", distance, "\"' will be deprecated; please use 'distance = \"gam\", link = \"", link, "\"' in the future."), call. = FALSE, immediate. = TRUE)
      distance <- "gam"
      attr(distance, "link") <- link
    }
    else if (tolower(distance) == "logit") {
      distance <- "glm"
      attr(distance, "link") <- "logit"
    }
    else if (!tolower(distance) %in% allowable.distances) {
      stop("The argument supplied to distance is not an allowable value. See ?matchit for allowable options.", call. = FALSE)
    }
    else {
      distance <- tolower(distance)
    }

  }
  else if (!is.numeric(distance) || (!is.null(dim(distance)) && length(dim(distance)) != 2)) {
    stop("'distance' must be a string with the name of the distance measure to be used or a numeric vector or matrix containing distance measures.", call. = FALSE)
  }
  else if (is.matrix(distance) && (is.null(method) || !method %in% c("nearest", "optimal", "full"))) {
    if (is.null(method)) method <- "NULL" else method <- paste0('"', method, '"')
    stop(paste0("'distance' cannot be supplied as a matrix with method = ", method, "."), call. = FALSE)
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
        stop("When supplied as a matrix, 'distance' must have dimensions NxN or N1xN0. See ?distance for details.", call. = FALSE)
      }
    }
    else {
      if (length(distance) != length(treat)) stop("'distance' must be the same length as the dataset if specified as a numeric vector.", call. = FALSE)
    }

    if (anyNA(distance)) stop("Missing values are not allowed in 'distance'.", call. = FALSE)
  }
  return(distance)
}

#Function to check ratio is acceptable
.process.ratio <- function(ratio = NULL, min.controls = NULL, max.controls = NULL, na.ok = FALSE) {
  if (length(ratio) == 0) ratio <- 1
  ratio.na <- anyNA(ratio)
  if (!na.ok && ratio.na) stop("'ratio' cannot be NA.", call. = FALSE)
  if (!ratio.na && (!is.atomic(ratio) || !is.numeric(ratio) || length(ratio) > 1 || ratio < 1)) {
    stop("'ratio' must be a single positive number.", call. = FALSE)
  }
  if (is.null(max.controls)) {
    ratio <- round(ratio)
    return(c(ratio = ratio))
  }
  else if (anyNA(max.controls) || !is.atomic(max.controls) || !is.numeric(max.controls) || length(max.controls) > 1) {
    stop("'max.controls' must be a single positive number.", call. = FALSE)
  }
  else if (ratio.na) {

  }
  else {
    if (!ratio.na && ratio <= 1) stop("'ratio' must be greater than 1 for variable ratio matching.", call. = FALSE)

    max.controls <- ceiling(max.controls)
    if (max.controls <= ratio) stop("'max.controls' must be greater than 'ratio' for variable ratio matching.", call. = FALSE)

    if (is.null(min.controls)) min.controls <- 1
    else if (anyNA(max.controls) || !is.atomic(max.controls) || !is.numeric(max.controls) || length(max.controls) > 1) {
        stop("'max.controls' must be a single positive number.", call. = FALSE)
    }
    else min.controls <- floor(min.controls)

    if (min.controls < 1) stop("'min.controls' cannot be less than 1 for variable ratio matching.", call. = FALSE)
    else if (min.controls >= ratio) stop("'min.controls' must be less than 'ratio' for variable ratio matching.", call. = FALSE)

    return(c(ratio = ratio, min.controls = min.controls, max.controls = max.controls))
  }
}
process.ratio <- function(ratio = NULL, method, min.controls = NULL, max.controls = NULL, ...) {
  #Should be run after process.inputs() and ignored inputs set to NULL
  ratio.null <- length(ratio) == 0
  ratio.na <- !ratio.null && anyNA(ratio)

  if (method %in% c("nearest", "optimal")) {
    if (ratio.null) ratio <- 1
    else if (ratio.na) stop("'ratio' cannot be NA.", call. = FALSE)
    else if (!is.atomic(ratio) || !is.numeric(ratio) || length(ratio) > 1 || ratio < 1) {
      stop("'ratio' must be a single positive number greater than or equal to 1.", call. = FALSE)
    }

    if (is.null(max.controls)) {
      ratio <- round(ratio)
    }
    else if (anyNA(max.controls) || !is.atomic(max.controls) || !is.numeric(max.controls) || length(max.controls) > 1) {
      stop("'max.controls' must be a single positive number.", call. = FALSE)
    }
    else {
      if (ratio <= 1) stop("'ratio' must be greater than 1 for variable ratio matching.", call. = FALSE)

      max.controls <- ceiling(max.controls)
      if (max.controls <= ratio) stop("'max.controls' must be greater than 'ratio' for variable ratio matching.", call. = FALSE)

      if (is.null(min.controls)) min.controls <- 1
      else if (anyNA(max.controls) || !is.atomic(max.controls) || !is.numeric(max.controls) || length(max.controls) > 1) {
        stop("'max.controls' must be a single positive number.", call. = FALSE)
      }
      else min.controls <- floor(min.controls)

      if (min.controls < 1) stop("'min.controls' cannot be less than 1 for variable ratio matching.", call. = FALSE)
      else if (min.controls >= ratio) stop("'min.controls' must be less than 'ratio' for variable ratio matching.", call. = FALSE)
    }
  }
  else if (method == "full") {
    if (!ratio.null) stop("'ratio' cannot be specified with full matching.", call. = FALSE) #should never be called because of process.inputs()
    if (!is.null(max.controls) && (anyNA(max.controls) || !is.atomic(max.controls) || !is.numeric(max.controls) || length(max.controls) > 1)) {
      stop("'max.controls' must be a single positive number.", call. = FALSE)
    }

    if (!is.null(min.controls) && (anyNA(min.controls) || !is.atomic(min.controls) || !is.numeric(min.controls) || length(min.controls) > 1)) {
      stop("'min.controls' must be a single positive number.", call. = FALSE)
    }
  }
  else if (method == "genetic") {
    if (ratio.null) ratio <- 1
    else if (ratio.na) stop("'ratio' cannot be NA.", call. = FALSE)
    else if (!is.atomic(ratio) || !is.numeric(ratio) || length(ratio) > 1 || ratio < 1) {
      stop("'ratio' must be a single positive number greater than or equal to 1.", call. = FALSE)
    }
    ratio <- round(ratio)

    min.controls <- max.controls <- NULL
  }
  else if (method == "cardinality") {
    if (ratio.null) ratio <- 1
    else if (!ratio.na && (!is.atomic(ratio) || !is.numeric(ratio) || length(ratio) > 1 || ratio < 0)) {
      stop("'ratio' must be a single positive number or NA.", call. = FALSE)
    }

    min.controls <- max.controls <- NULL
  }

  return(list(ratio = ratio, min.controls = min.controls, max.controls = max.controls))

}

#Function to check if caliper is okay and process it
process.caliper <- function(caliper = NULL, method, data = NULL, covs = NULL, mahcovs = NULL, distance = NULL, discarded = NULL, std.caliper = TRUE) {

  #Check method; must be able to use a caliper
  #Check caliper names; if "" is one of them but distance = "mahal", throw error;
  #otherwise make sure variables exist in data or covs
  #Make sure no calipers are used on binary or factor variables (throw error if so)
  #Ignore calipers used on single-value variables or with caliper = NA or Inf
  #Export caliper.formula to add to covs
  #If std, export standardized versions

  #Check need for caliper
  if (length(caliper) == 0 || is.null(method) || !method %in% c("nearest", "genetic", "full")) return(NULL)

  #Check if form of caliper is okay
  if (!is.atomic(caliper) || !is.numeric(caliper)) stop("'caliper' must be a numeric vector.", call. = FALSE)

  #Check caliper names
  if (length(caliper) == 1 && (is.null(names(caliper)) || identical(names(caliper), ""))) names(caliper) <- ""
  else if (is.null(names(caliper))) stop("'caliper' must be a named vector with names corresponding to the variables for which a caliper is to be applied.", call. = FALSE)
  else if (anyNA(names(caliper))) stop("'caliper' names cannot include NA.", call. = FALSE)
  else if (sum(names(caliper) == "") > 1) stop("No more than one entry in 'caliper' can have no name.", call. = FALSE)

  if (any(names(caliper) == "") && is.null(distance)) stop("All entries in 'caliper' must be named when distance = \"mahalanobis\".", call. = FALSE)

  #Check if caliper name is in available data
  cal.in.data <- setNames(names(caliper) %in% names(data), names(caliper))
  cal.in.covs <- setNames(names(caliper) %in% names(covs), names(caliper))
  cal.in.mahcovs <- setNames(names(caliper) %in% names(mahcovs), names(caliper))
  if (any(names(caliper) != "" & !cal.in.covs & !cal.in.data)) stop(paste0("All variables named in 'caliper' must be in 'data'. Variables not in 'data':\n\t",
                                         paste0(names(caliper)[names(caliper) != "" & !cal.in.data & !cal.in.covs & !cal.in.mahcovs], collapse = ", ")), call. = FALSE)

  #Check std.caliper
  if (length(std.caliper) == 0 || !is.atomic(std.caliper) || !is.logical(std.caliper)) stop("'std.caliper' must be a logical (TRUE/FALSE) vector.", call. = FALSE)
  if (length(std.caliper) == 1) std.caliper <- setNames(rep.int(std.caliper, length(caliper)), names(caliper))
  else if (length(std.caliper) != length(caliper)) stop("'std.caliper' must be the same length as 'caliper'", call. = FALSE)
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
    else {
      if (x == "") var <- distance
      else if (cal.in.data[x]) var <- data[[x]]
      else if (cal.in.covs[x]) var <- covs[[x]]
      else var <- mahcovs[[x]]

      return(is.factor(var) || is.character(var))
    }
  }, logical(1L))

  if (any(cat.vars)) {
    stop(paste0("Calipers cannot be used with binary, factor, or character variables. Offending variables:\n\t",
                paste0(ifelse(names(caliper) == "", "<distance>", names(caliper))[cat.vars], collapse = ", ")), call. = FALSE)
  }

  #Process calipers according to std.caliper
  std.caliper <- std.caliper[names(std.caliper) %in% names(caliper)]
  if (anyNA(std.caliper)) stop("'std.caliper' cannot be NA.", call. = FALSE)

  if (any(std.caliper)) {
    if ("" %in% names(std.caliper) && isTRUE(std.caliper[names(std.caliper) == ""]) && is.matrix(distance)) {
      stop("When 'distance' is supplied as a matrix and a caliper for it is specified, 'std.caliper' must be FALSE for the distance measure.", call. = FALSE)
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

  return(abs(caliper))

}

#Function to ensure no subclass is devoid of both treated and control units by "scooting" units
#from other subclasses. From WeightIt.
subclass_scoot <- function(sub, treat, x, min.n = 1) {
  #Reassigns subclasses so there are no empty subclasses
  #for each treatment group. Copied from WeightIt with
  #slight modifications.

  treat <- as.character(treat)
  unique.treat <- unique(treat, nmax = 2)

  names(x) <- seq_along(x)
  names(sub) <- seq_along(sub)
  original.order <- names(x)

  nsub <- length(unique(sub))

  #Turn subs into a contiguous sequence
  sub <- setNames(setNames(seq_len(nsub), sort(unique(sub)))[as.character(sub)],
                  original.order)

  if (any(table(treat) < nsub * min.n)) {
    stop(paste0("Not enough units to fit ", min.n, ngettext(min.n, " treated and control unit",
                                                            " treated and control units"),
                " in each subclass."), call. = FALSE)
  }

  for (t in unique.treat) {
    if (length(x[treat == t]) == nsub) {
      sub[treat == t] <- seq_len(nsub)
    }
  }

  sub_tab <- table(treat, sub)

  if (any(sub_tab < min.n)) {

    soft_thresh <- function(x, minus = 1) {
      x <- x - minus
      x[x < 0] <- 0
      x
    }

    for (t in unique.treat) {
      for (n in seq_len(min.n)) {
        while (any(sub_tab[t,] == 0)) {
          first_0 <- which(sub_tab[t,] == 0)[1]

          if (first_0 == nsub ||
              (first_0 != 1 &&
               sum(soft_thresh(sub_tab[t, seq(1, first_0 - 1)]) / abs(first_0 - seq(1, first_0 - 1))) >=
               sum(soft_thresh(sub_tab[t, seq(first_0 + 1, nsub)]) / abs(first_0 - seq(first_0 + 1, nsub))))) {
            #If there are more and closer nonzero subs to the left...
            first_non0_to_left <- max(seq(1, first_0 - 1)[sub_tab[t, seq(1, first_0 - 1)] > 0])

            name_to_move <- names(sub)[which(x == max(x[treat == t & sub == first_non0_to_left]) & treat == t & sub == first_non0_to_left)[1]]

            sub[name_to_move] <- first_0
            sub_tab[t, first_0] <- 1L
            sub_tab[t, first_non0_to_left] <- sub_tab[t, first_non0_to_left] - 1L

          }
          else {
            #If there are more and closer nonzero subs to the right...
            first_non0_to_right <- min(seq(first_0 + 1, nsub)[sub_tab[t, seq(first_0 + 1, nsub)] > 0])

            name_to_move <- names(sub)[which(x == min(x[treat == t & sub == first_non0_to_right]) & treat == t & sub == first_non0_to_right)[1]]

            sub[name_to_move] <- first_0
            sub_tab[t, first_0] <- 1L
            sub_tab[t, first_non0_to_right] <- sub_tab[t, first_non0_to_right] - 1L
          }
        }

        sub_tab[t,] <- sub_tab[t,] - 1
      }
    }

    #Unsort
    sub <- sub[names(sub)]
  }

  return(sub)
}

#Function to check if package is installed. From WeightIt.
check.package <- function(package.name, alternative = FALSE) {
  packages.not.installed <- package.name[!vapply(package.name, requireNamespace, logical(1L),
                                                 quietly = TRUE)]
  if (length(packages.not.installed) > 0) {
    if (alternative) return(FALSE)
    else {
      plural <- length(packages.not.installed) > 1
      stop(paste0("Package", if (plural) "s " else " ",
                  word_list(packages.not.installed, quotes = 1, is.are = TRUE),
                  " needed for this function to work. Please install ",
                  if (plural) "them" else "it","."),
           call. = FALSE)
    }
  }
  else return(invisible(TRUE))
}

#Create info component of matchit object
create_info <- function(method, fn1, link, discard, replace, ratio, max.controls, mcall, mahalanobis, subclass, antiexact, distance_is_matrix) {
  info <- list(method = method,
               distance = if (is.null(fn1)) NULL else sub("distance2", "", fn1, fixed = TRUE),
               link = if (is.null(link)) NULL else link,
               discard = discard,
               replace = if (!is.null(method) && method %in% c("nearest", "genetic")) replace else NULL,
               ratio = if (!is.null(method) && method %in% c("nearest", "optimal", "genetic")) ratio else NULL,
               max.controls = if (!is.null(method) && method %in% c("nearest", "optimal")) max.controls else NULL,
               # mahalanobis = is.full.mahalanobis || !is.null(mahvars),
               mahalanobis = mahalanobis,
               # subclass = if (!is.null(method) && method == "subclass") length(unique(match.out$subclass[!is.na(match.out$subclass)])) else NULL,
               subclass = if (!is.null(method) && method == "subclass") length(unique(subclass[!is.na(subclass)])) else NULL,
               # antiexact = if (!is.null(antiexact)) colnames(antiexactcovs) else NULL
               antiexact = antiexact,
               distance_is_matrix = distance_is_matrix)
}

#Function to turn a method name into a phrase describing the method
info.to.method <- function(info) {

  out.list <- setNames(vector("list", 3), c("kto1", "type", "replace"))
  out.list[["kto1"]] <- if (!is.null(info$ratio)) paste0(if (!is.null(info$max.controls)) "variable ratio ", round(info$ratio, 2), ":1") else NULL
  out.list[["type"]] <- if (is.null(info$method)) "none (no matching)" else
    switch(info$method,
           "exact" = "exact matching",
           "cem" = "coarsened exact matching",
           "nearest" = "nearest neighbor matching",
           "optimal" = "optimal pair matching",
           "full" = "optimal full matching",
           "genetic" = "genetic matching",
           "subclass" = paste0("subclassification (", info$subclass, " subclasses)"),
           "cardinality" = "cardinality matching")
  out.list[["replace"]] <- if (!is.null(info$replace) && info$method %in% c("nearest", "optimal", "genetic")) {
    if (info$replace) "with replacement"
    else "without replacement"
  } else NULL

  firstup(do.call("paste", c(unname(out.list), list(sep = " "))))
}

info.to.distance <- function(info) {
  distance <- info$distance
  link <- info$link
  if (!is.null(link) && startsWith(as.character(link), "linear")) {
    linear <- TRUE
    link <- sub("linear.", "", as.character(link))
  }
  else linear <- FALSE

  if (distance == "glm") {
    if (link == "logit") dist <- "logistic regression"
    else if (link == "probit") dist <- "probit regression"
    else dist <- paste("GLM with a", link, "link")
  }
  else if (distance == "gam") {
    dist <- paste("GAM with a", link, "link")
  }
  else if (distance == "rpart") {
    dist <- "CART"
  }
  else if (distance == "nnet") {
    dist <- "a neural network"
  }
  else if (distance == "cbps") {
    dist <- "CBPS"
  }
  else if (distance == "bart") {
    dist <- "BART"
  }
  else if (distance == "randomforest") {
    dist <- "a random forest"
  }

  if (linear) dist <- paste(dist, "and linearized")

  return(dist)
}

#Function to turn a vector into a string with "," and "and" or "or" for clean messages. 'and.or'
#controls whether words are separated by "and" or "or"; 'is.are' controls whether the list is
#followed by "is" or "are" (to avoid manually figuring out if plural); quotes controls whether
#quotes should be placed around words in string. From WeightIt.
word_list <- function(word.list = NULL, and.or = c("and", "or"), is.are = FALSE, quotes = FALSE) {
  #When given a vector of strings, creates a string of the form "a and b"
  #or "a, b, and c"
  #If is.are, adds "is" or "are" appropriately
  L <- length(word.list)
  word.list <- add_quotes(word.list, quotes)

  if (L == 0) {
    out <- ""
    attr(out, "plural") = FALSE
  }
  else {
    word.list <- word.list[!word.list %in% c(NA_character_, "")]
    L <- length(word.list)
    if (L == 0) {
      out <- ""
      attr(out, "plural") = FALSE
    }
    else if (L == 1) {
      out <- word.list
      if (is.are) out <- paste(out, "is")
      attr(out, "plural") = FALSE
    }
    else {
      and.or <- match_arg(and.or)
      if (L == 2) {
        out <- paste(word.list, collapse = paste0(" ", and.or," "))
      }
      else {
        out <- paste(paste(word.list[seq_len(L-1)], collapse = ", "),
                     word.list[L], sep = paste0(", ", and.or," "))

      }
      if (is.are) out <- paste(out, "are")
      attr(out, "plural") = TRUE
    }

  }
  return(out)
}

#Add quotation marks around a string.
add_quotes <- function(x, quotes = 2) {
  if (!isFALSE(quotes)) {
    if (isTRUE(quotes) || as.integer(quotes) == 2) x <- paste0("\"", x, "\"")
    else if (as.integer(quotes) == 1) x <- paste0("\'", x, "\'")
    else stop("'quotes' must be boolean, 1, or 2.")
  }
  x
}

#More informative and cleaner version of base::match.arg. From WeightIt.
match_arg <- function(arg, choices, several.ok = FALSE) {
  #Replaces match.arg() but gives cleaner error message and processing
  #of arg.
  if (missing(arg))
    stop("No argument was supplied to match_arg.", call. = FALSE)
  arg.name <- paste(deparse(substitute(arg), width.cutoff = 500L), collapse = " ")

  if (missing(choices)) {
    formal.args <- formals(sys.function(sysP <- sys.parent()))
    choices <- eval(formal.args[[as.character(substitute(arg))]],
                    envir = sys.frame(sysP))
  }

  if (is.null(arg))
    return(choices[1L])
  else if (!is.character(arg))
    stop(paste0("The argument to '", arg.name, "' must be NULL or a character vector"), call. = FALSE)
  if (!several.ok) {
    if (identical(arg, choices))
      return(arg[1L])
    if (length(arg) > 1L)
      stop(paste0("The argument to '", arg.name, "' must be of length 1"), call. = FALSE)
  }
  else if (length(arg) == 0)
    stop(paste0("The argument to '", arg.name, "' must be of length >= 1"), call. = FALSE)

  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
  if (all(i == 0L))
    stop(paste0("The argument to '", arg.name, "' should be ", if (length(choices) > 1) {if (several.ok) "at least one of " else "one of "} else "",
                word_list(choices, and.or = "or", quotes = 2), "."),
         call. = FALSE)
  i <- i[i > 0L]
  if (!several.ok && length(i) > 1)
    stop("There is more than one match in 'match_arg'")
  choices[i]
}

#Turn a vector into a 0/1 vector. 'zero' and 'one' can be supplied to make it clear which is
#which; otherwise, a guess is used. From WeightIt.
binarize <- function(variable, zero = NULL, one = NULL) {
  if (length(unique(variable)) > 2) stop(paste0("Cannot binarize ", paste(deparse(substitute(variable)), collapse = " "), ": more than two levels."))
  if (is.character(variable) || is.factor(variable)) {
    variable <- factor(variable, nmax = 2)
    unique.vals <- levels(variable)
  }
  else {
    unique.vals <- unique(variable, nmax = 2)
  }

  if (is.null(zero)) {
    if (is.null(one)) {
      if (can_str2num(unique.vals)) {
        variable.numeric <- str2num(variable)
      }
      else {
        variable.numeric <- as.numeric(variable)
      }

      if (0 %in% variable.numeric) zero <- 0
      else zero <- min(variable.numeric, na.rm = TRUE)

      return(setNames(as.integer(variable.numeric != zero), names(variable)))
    }
    else {
      if (one %in% unique.vals) return(setNames(as.integer(variable == one), names(variable)))
      else stop("The argument to 'one' is not the name of a level of variable.", call. = FALSE)
    }
  }
  else {
    if (zero %in% unique.vals) return(setNames(as.integer(variable != zero), names(variable)))
    else stop("The argument to 'zero' is not the name of a level of variable.", call. = FALSE)
  }
}

#Make interaction vector out of matrix of covs
exactify <- function(X, nam = NULL, sep = "|", include_vars = FALSE) {
  if (is.null(nam)) nam <- rownames(X)
  if (is.matrix(X)) X <- setNames(lapply(seq_len(ncol(X)), function(i) X[,i]), colnames(X))
  if (!is.list(X)) stop("X must be a matrix, data frame, or list.")

  #Ensure no ambiguity is created by sep
  sep0 <- sep
  unique.x <- unlist(lapply(X, function(x) as.character(unique(x))))
  while (any(grepl(sep, unique.x, fixed = TRUE))) {
    sep0 <- paste0(sep0, sep)
  }

  if (include_vars) {
    for (i in seq_along(X)) {
      if (is.character(X[[i]]) || is.factor(X[[i]])) {
        X[[i]] <- paste0(names(X)[i], ' = "', X[[i]], '"')
      }
      else {
        X[[i]] <- paste0(names(X)[i], ' = ', X[[i]])
      }
    }
  }

  out <- do.call("paste", c(X, sep = sep0))
  if (!is.null(nam)) names(out) <- nam
  out
}

#Determine whether a character vector can be coerced to numeric
can_str2num <- function(x) {
  nas <- is.na(x)
  suppressWarnings(x_num <- as.numeric(as.character(x[!nas])))
  return(!anyNA(x_num))
}

#Cleanly coerces a character vector to numeric; best to use after can_str2num()
str2num <- function(x) {
  nas <- is.na(x)
  suppressWarnings(x_num <- as.numeric(as.character(x)))
  x_num[nas] <- NA
  return(x_num)
}

#Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

#Clean printing of data frames with numeric and NA elements.
round_df_char <- function(df, digits, pad = "0", na_vals = "") {
  #Digits is passed to round(). pad is used to replace trailing zeros so decimal
  #lines up. Should be "0" or " "; "" (the empty string) un-aligns decimals.
  #na_vals is what NA should print as.

  if (NROW(df) == 0 || NCOL(df) == 0) return(df)
  if (!is.data.frame(df)) df <- as.data.frame.matrix(df, stringsAsFactors = FALSE)
  rn <- rownames(df)
  cn <- colnames(df)

  infs <- o.negs <- array(FALSE, dim = dim(df))
  nas <- is.na(df)
  nums <- vapply(df, is.numeric, logical(1))
  infs[,nums] <- vapply(which(nums), function(i) !nas[,i] & !is.finite(df[[i]]), logical(NROW(df)))

  for (i in which(!nums)) {
    if (can_str2num(df[[i]])) {
      df[[i]] <- str2num(df[[i]])
      nums[i] <- TRUE
    }
  }

  o.negs[,nums] <- !nas[,nums] & df[nums] < 0 & round(df[nums], digits) == 0
  df[nums] <- round(df[nums], digits = digits)

  for (i in which(nums)) {
    df[[i]] <- format(df[[i]], scientific = FALSE, justify = "none", trim = TRUE,
                      drop0trailing = !identical(as.character(pad), "0"))

    if (!identical(as.character(pad), "0") && any(grepl(".", df[[i]], fixed = TRUE))) {
      s <- strsplit(df[[i]], ".", fixed = TRUE)
      lengths <- lengths(s)
      digits.r.of.. <- rep(0, NROW(df))
      digits.r.of..[lengths > 1] <- nchar(vapply(s[lengths > 1], `[[`, character(1L), 2))
      max.dig <- max(digits.r.of..)

      dots <- ifelse(lengths > 1, "", if (as.character(pad) != "") "." else pad)
      pads <- vapply(max.dig - digits.r.of.., function(n) paste(rep(pad, n), collapse = ""), character(1L))

      df[[i]] <- paste0(df[[i]], dots, pads)
    }
  }

  df[o.negs] <- paste0("-", df[o.negs])

  # Insert NA placeholders
  df[nas] <- na_vals
  df[infs] <- "N/A"

  if (length(rn) > 0) rownames(df) <- rn
  if (length(cn) > 0) names(df) <- cn

  return(df)
}

#Generalized inverse; port of MASS::ginv()
generalized_inverse <-function(sigma) {
  sigmasvd <- svd(sigma)
  pos <- sigmasvd$d > max(1e-8 * sigmasvd$d[1L], 0)
  sigma_inv <- sigmasvd$v[, pos, drop = FALSE] %*% (sigmasvd$d[pos]^-1 * t(sigmasvd$u[, pos, drop = FALSE]))
  return(sigma_inv)
}

#Get covariates (RHS) vars from formula
get.covs.matrix <- function(formula = NULL, data = NULL) {

  if (is.null(formula)) {
    fnames <- colnames(data)
    fnames[!startsWith(fnames, "`")] <- paste0("`", fnames[!startsWith(fnames, "`")], "`")
    formula <- reformulate(fnames)
  }
  else formula <- update(terms(formula, data = data), NULL ~ . + 1)

  mf <- model.frame(terms(formula, data = data), data,
                    na.action = na.pass)

  chars.in.mf <- vapply(mf, is.character, logical(1L))
  mf[chars.in.mf] <- lapply(mf[chars.in.mf], factor)

  X <- model.matrix(formula, data = mf,
                    contrasts.arg = lapply(Filter(is.factor, mf),
                                           contrasts, contrasts = FALSE))
  assign <- attr(X, "assign")[-1]
  X <- X[,-1,drop=FALSE]
  attr(X, "assign") <- assign

  return(X)
}

#Extracts and names the "assign" attribute from get.covs.matrix()
get_assign <- function(mat) {
  if (is.null(attr(mat, "assign"))) return(NULL)

  setNames(attr(mat, "assign"), colnames(mat))
}

#Convert match.matrix (mm) using numerical indices to using char rownames
nummm2charmm <- function(nummm, treat) {
  #Assumes nummm has rownames
  charmm <- matrix(NA_character_, nrow = nrow(nummm), ncol = ncol(nummm),
                   dimnames = dimnames(nummm))
  charmm[] <- names(treat)[nummm]
  charmm
}

charmm2nummm <- function(charmm, treat) {
  nummm <- matrix(NA_integer_, nrow = nrow(charmm), ncol = ncol(charmm))
  n_index <- setNames(seq_along(treat), names(treat))
  nummm[] <- n_index[charmm]
  nummm
}

#Get subclass from match.matrix. Only to be used if replace = FALSE. See subclass2mmC.cpp for reverse.
mm2subclass <- function(mm, treat) {
  lab <- names(treat)
  ind1 <- which(treat == 1)

  subclass <- setNames(rep(NA_character_, length(treat)), lab)
  no.match <- is.na(mm)
  subclass[ind1[!no.match[,1]]] <- ind1[!no.match[,1]]
  subclass[mm[!no.match]] <- ind1[row(mm)[!no.match]]

  subclass <- setNames(factor(subclass, nmax = length(ind1)), lab)
  levels(subclass) <- seq_len(nlevels(subclass))

  return(subclass)
}

#(Weighted) variance that uses special formula for binary variables
wvar <- function(x, bin.var = NULL, w = NULL) {
  if (is.null(w)) w <- rep(1, length(x))
  if (is.null(bin.var)) bin.var <- all(x == 0 | x == 1)

  w <- w / sum(w) #weights normalized to sum to 1
  mx <- sum(w * x) #weighted mean

  if (bin.var) {
    mx*(1-mx)
  }
  else {
    #Reliability weights variance; same as cov.wt()
    sum(w * (x - mx)^2)/(1 - sum(w^2))
  }
}

#Weighted mean faster than weighted.mean()
wm <- function(x, w = NULL, na.rm = TRUE) {
  if (is.null(w)) {
    if (anyNA(x)) {
      if (!na.rm) return(NA_real_)
      nas <- which(is.na(x))
      x <- x[-nas]
    }
    return(sum(x)/length(x))
  }
  else {
    if (anyNA(x) || anyNA(w)) {
      if (!na.rm) return(NA_real_)
      nas <- which(is.na(x) | is.na(w))
      x <- x[-nas]
      w <- w[-nas]
    }
    return(sum(x*w)/sum(w))
  }
}

#Effective sample size
ESS <- function(w) {
  sum(w)^2/sum(w^2)
}

#Compute sample sizes
nn <- function(treat, weights, discarded, s.weights) {
  if (is.null(s.weights)) s.weights <- rep(1, length(treat))
  weights <- weights * s.weights
  n <- matrix(0, ncol=2, nrow=6, dimnames = list(c("All (ESS)", "All", "Matched (ESS)","Matched", "Unmatched","Discarded"),
                                                  c("Control", "Treated")))

  #                       Control                                    Treated
  n["All (ESS)",] <-     c(ESS(s.weights[treat==0]),                ESS(s.weights[treat==1]))
  n["All",] <-           c(sum(treat==0),                           sum(treat==1))
  n["Matched (ESS)",] <- c(ESS(weights[treat==0]),                  ESS(weights[treat==1]))
  n["Matched",] <-       c(sum(treat==0 & weights > 0),             sum(treat==1 & weights > 0))
  n["Unmatched",] <-     c(sum(treat==0 & weights==0 & !discarded), sum(treat==1 & weights==0 & !discarded))
  n["Discarded",] <-     c(sum(treat==0 & discarded),               sum(treat==1 & discarded))

  return(n)
}

#Compute subclass sample sizes
qn <- function(treat, subclass, discarded) {

  qn <- table(treat[!discarded], subclass[!discarded])
  dimnames(qn) <- list(c("Control", "Treated"), levels(subclass))

  if (any(discarded)) {
    qn <- cbind(qn, table(treat[discarded]))
    colnames(qn)[ncol(qn)] <- "Discarded"
  }
  qn <- rbind(qn, colSums(qn))
  rownames(qn)[nrow(qn)] <- "Total"

  qn <- cbind(qn, rowSums(qn))
  colnames(qn)[ncol(qn)] <- "All"

  return(qn)
}

#Used to load backports functions. No need to touch, but must always be included somewhere.
.onLoad <- function(libname, pkgname) {
  backports::import(pkgname)
}