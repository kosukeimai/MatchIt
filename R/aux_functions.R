check.inputs <- function(method, distance, mcall, exact, mahvars, caliper, discard, replace, ratio, m.order) {
  method <- match_arg(method, c("exact", "cem", "nearest", "optimal", "full", "genetic", "subclass"))

  ignored.inputs <- character(0)
  error.inputs <- character(0)
  if (method == "exact") {
    for (i in c("distance", "exact", "mahvars", "caliper", "discard", "replace", "ratio", "m.order")) {
      if (i %in% names(mcall) && !is.null(get0(i))) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "cem") {
    for (i in c("distance", "exact", "mahvars", "caliper", "discard", "replace", "ratio", "m.order")) {
      if (i %in% names(mcall) && !is.null(get0(i))) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "nearest") {
    if (distance == "mahalanobis") {
      for (e in c("mahvars", "caliper", "discard")) {
        if (e %in% names(mcall) && !is.null(get0(e))) {
          error.inputs <- c(error.inputs, e)
        }
      }
    }
    else {
      if (!is.null(mahvars) && is.null(caliper)) {
        warning("When mahvars are specified, an argument should be supplied to caliper.", call. = FALSE)
      }
    }
  }
  else if (method == "optimal") {
    if (distance == "mahalanobis") {
      for (e in c("mahvars", "caliper", "discard")) {
        if (e %in% names(mcall) && !is.null(get0(e))) {
          error.inputs <- c(error.inputs, e)
        }
      }
    }
    else {
      if (!is.null(mahvars) && is.null(caliper)) {
        warning("When mahvars are specified, an argument should be supplied to caliper.", call. = FALSE)
      }
    }

    for (i in c("replace", "m.order")) {
      if (i %in% names(mcall) && !is.null(get0(i))) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }

  }
  else if (method == "full") {
    if (distance == "mahalanobis") {
      for (e in c("mahvars", "caliper", "discard")) {
        if (e %in% names(mcall) && !is.null(get0(e))) {
          error.inputs <- c(error.inputs, e)
        }
      }
    }
    else {
      if (!is.null(mahvars) && is.null(caliper)) {
        warning("When mahvars are specified, an argument should be supplied to caliper.", call. = FALSE)
      }
    }

    for (i in c("replace", "ratio", "m.order")) {
      if (i %in% names(mcall) && !is.null(get0(i))) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }
  else if (method == "genetic") {
    if (distance == "mahalanobis") {
      for (e in c("mahvars", "caliper", "discard")) {
        if (e %in% names(mcall) && !is.null(get0(e))) {
          error.inputs <- c(error.inputs, e)
        }
      }
    }
    else {
      if (!is.null(mahvars) && is.null(caliper)) {
        warning("When mahvars are specified, an argument should be supplied to caliper.", call. = FALSE)
      }
    }
  }
  else if (method == "subclass") {
    if (distance == "mahalanobis") {
      stop("distance = \"mahalanobis\" is not compatible with subclassification.", call. = FALSE)
    }

    for (i in c("exact", "mahvars", "caliper", "replace", "ratio", "m.order")) {
      if (i %in% names(mcall) && !is.null(get0(i))) {
        ignored.inputs <- c(ignored.inputs, i)
      }
    }
  }

  if (length(ignored.inputs) > 0) warning(paste0(ngettext(length(ignored.inputs), "The argument ", "The arguments "),
                                                 word_list(ignored.inputs, quotes = 1, is.are = TRUE),
                                                 " not used with method = \"", method, "\" and will be ignored."),
                                          call. = FALSE, immediate. = TRUE)
  if (length(error.inputs) > 0) stop(paste0(ngettext(length(error.inputs), "The argument ", "The arguments "),
                                            word_list(error.inputs, quotes = 1, is.are = TRUE),
                                            " not allowed with method = \"", method,
                                            "\" and distance = \"", distance, "\"."),
                                     call. = FALSE, immediate. = TRUE)
}

process.distance <- function(distance, method) {
  if (is.null(distance)) stop(paste0("distance cannot be NULL with method = \"", method, "\"."), call. = FALSE)
  else if (is.vector(distance, "character") && length(distance) == 1) {
    allowable.distances <- c("glm", "cbps", "gam", "mahalanobis", "nnet", "rpart", "bart")

    if (tolower(distance) %in% c("cauchit", "cloglog", "linear.cloglog", "linear.log", "linear.logit", "linear.probit",
                        "linear.cauchit", "log", "probit")) {
      warning(paste0("'distance = \"", distance, "\"' will be deprecated; please use 'distance = \"glm\", link = \"", distance, "\"' in the future."), call. = FALSE, immediate. = TRUE)
      link <- distance
      distance <- "glm"
      attr(distance, "link") <- link
    }
    else if (tolower(distance) %in% tolower(c("GAMcloglog", "GAMlog", "GAMlogit", "GAMprobit"))) {
      warning(paste0("'distance = \"", distance, "\"' will be deprecated; please use 'distance = \"gam\", link = \"", sub("GAM", "", distance), "\"' in the future."), call. = FALSE, immediate. = TRUE)
      link <- distance
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
  else if (!is.vector(distance, "numeric")) {
    stop("distance must be a string with the name of the distance measure to be used or a numeric vector containing distance measures.", call. = FALSE)
  }
  return(distance)
}

process.ratio <- function(ratio) {
  if (length(ratio) == 0) ratio <- 1
  if (!is.atomic(ratio) || !is.numeric(ratio) || length(ratio) > 1 || ratio < 1) {
    stop("Ratio must be a single positive number.", call. = FALSE)
  }
  round(ratio)
}

subclass_scoot <- function(sub, treat, x) {
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

  if (any(table(treat) < nsub)) {
    stop("Too many subclasses were requested.", call. = FALSE)
  }

  for (t in unique.treat) {
    if (length(x[treat == t]) == nsub) {
      sub[treat == t] <- seq_len(nsub)
    }
  }

  sub_tab <- table(treat, sub)

  if (any(sub_tab == 0)) {

    soft_thresh <- function(x, minus = 1) {
      x <- x - minus
      x[x < 0] <- 0
      x
    }

    for (t in unique.treat) {
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
    }

    #Unsort
    sub <- sub[names(sub)]
  }

  return(sub)
}

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

word_list <- function(word.list = NULL, and.or = c("and", "or"), is.are = FALSE, quotes = FALSE) {
  #When given a vector of strings, creates a string of the form "a and b"
  #or "a, b, and c"
  #If is.are, adds "is" or "are" appropriately
  L <- length(word.list)
  if (quotes) {
    if (as.integer(quotes) == 2) word.list <- vapply(word.list, function(x) paste0("\"", x, "\""), character(1L))
    else if (as.integer(quotes) == 1) word.list <- vapply(word.list, function(x) paste0("\'", x, "\'"), character(1L))
    else stop("'quotes' must be boolean, 1, or 2.")
  }
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

match_arg <- function(arg, choices, several.ok = FALSE) {
  #Replaces match.arg() but gives cleaner error message and processing
  #of arg.
  if (missing(arg))
    stop("No argument was supplied to match_arg.", call. = FALSE)
  arg.name <- deparse1(substitute(arg))

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

.onLoad <- function(libname, pkgname) {
  backports::import(pkgname)
}