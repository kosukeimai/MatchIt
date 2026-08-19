#Add quotes to a string
add_quotes <- function(x, quotes = 2L) {
  if (isFALSE(quotes)) {
    return(x)
  }

  if (isTRUE(quotes)) {
    quotes <- '"'
  }

  if (rlang::is_string(quotes)) {
    return(paste0(quotes, x, str_rev(quotes)))
  }

  if (!rlang::is_integerish(quotes) || quotes > 2L) {
    stop("`quotes` must be boolean, 1, 2, or a string.")
  }

  if (quotes == 0L) {
    return(x)
  }

  x <- {
    if (quotes == 1L) sprintf("'%s'", x)
    else sprintf('"%s"', x)
  }

  x
}

# Version of interaction(., drop = TRUE) that doesn't succumb to vector limit reached by
# avoiding Cartesian expansion. Falls back to interaction() for small problems.
interaction2 <- function(..., sep = ".", lex.order = TRUE) {

  narg <- ...length()

  if (narg == 0L) {
    stop("No factors specified")
  }

  if (narg == 1L && is.list(..1)) {
    args <- ..1
    narg <- length(args)
  }
  else {
    args <- list(...)
  }

  for (i in seq_len(narg)) {
    args[[i]] <- as.factor(args[[i]])
  }

  if (do.call("prod", lapply(args, nlevels)) <= 1e6) {
    return(interaction(args, drop = TRUE, sep = sep,
                       lex.order = if (is.null(lex.order)) TRUE else lex.order))
  }

  out <- do.call(function(...) paste(..., sep = sep), args)

  args_char <- lapply(args, function(x) {
    x <- unclass(x)
    formatC(x, format = "d", flag = "0", width = max(1, ceiling(log10(max(x)))))
  })

  lev <- {
    if (is.null(lex.order)) unique(out)
    else if (lex.order) unique(out[order(do.call("paste", c(args_char, list(sep = sep))))])
    else unique(out[order(do.call("paste", c(rev(args_char), list(sep = sep))))])
  }

  factor(out, levels = lev)
}

#Turn a vector into a 0/1 vector. 'zero' and 'one' can be supplied to make it clear which is
#which; otherwise, a guess is used. From WeightIt.
binarize <- function(variable, zero = NULL, one = NULL) {
  var.name <- deparse1(substitute(variable))
  if (is.character(variable) || is.factor(variable)) {
    variable <- factor(variable, nmax = if (is.factor(variable)) nlevels(variable) else NA)
    unique.vals <- levels(variable)
  }
  else {
    unique.vals <- unique(variable)
  }

  if (length(unique.vals) == 1L) {
    return(rep_with(1L, variable))
  }

  if (length(unique.vals) != 2L) {
    arg::err("cannot binarize {.var {var.name}}: more than two levels")
  }

  if (is_not_null(zero)) {
    if (!zero %in% unique.vals) {
      arg::err("the argument to {.arg zero} is not the name of a level of {.var {var.name}}")
    }

    return(setNames(as.integer(variable != zero), names(variable)))
  }

  if (is_not_null(one)) {
    if (!one %in% unique.vals) {
      arg::err("the argument to {.arg one} is not the name of a level of {.var {var.name}}")
    }

    return(setNames(as.integer(variable == one), names(variable)))
  }

  if (is.logical(variable)) {
    return(setNames(as.integer(variable), names(variable)))
  }

  if (is.numeric(variable)) {
    zero <- {
      if (any(unique.vals == 0)) 0
      else min(unique.vals, na.rm = TRUE)
    }

    return(setNames(as.integer(variable != zero), names(variable)))
  }

  variable.numeric <- {
    if (can_str2num(unique.vals)) setNames(str2num(unique.vals), unique.vals)[variable]
    else as.numeric(factor(variable, levels = unique.vals))
  }

  zero <- {
    if (0 %in% variable.numeric) 0
    else min(variable.numeric, na.rm = TRUE)
  }

  setNames(as.integer(variable.numeric != zero), names(variable))
}

is_null <- function(x) {isTRUE(length(x) == 0L)}
is_not_null <- function(x) !is_null(x)
`%or%` <- function(x, y) {
  # like `%||%` but works for non-NULL length 0 objects
  if (is_null(x)) y else x
}
allNA <- function(x) {
  anyNA(x) && all(is.na(x))
}

null_or_error <- function(x) {is_null(x) || inherits(x, "try-error")}

#Determine whether a character vector can be coerced to numeric
can_str2num <- function(x) {
  if (is.numeric(x) || is.logical(x)) {
    return(TRUE)
  }

  nas <- is.na(x)
  x_num <- suppressWarnings(as.numeric(as.character(x[!nas])))

  !anyNA(x_num)
}

#Cleanly coerces a character vector to numeric; best to use after can_str2num()
str2num <- function(x) {
  nas <- is.na(x)
  if (!is.numeric(x) && !is.logical(x)) x <- as.character(x)
  x_num <- suppressWarnings(as.numeric(x))
  is.na(x_num[nas]) <- TRUE
  x_num
}

#Capitalize first letter of string
firstup <- function(x) {
  substr(x, 1L, 1L) <- toupper(substr(x, 1L, 1L))
  x
}

#Capitalize first letter of each word
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste0(toupper(substring(s, 1L, 1L)),
                            {s <- substring(s, 2L)
                            if (strict) tolower(s) else s},
                            collapse = " ")
  sapply(strsplit(s, split = " ", fixed = TRUE), cap,
         USE.NAMES = is_not_null(names(s)))
}

#Reverse a string
str_rev <- function(x) {
  vapply(lapply(strsplit(x, NULL), rev), paste, character(1L), collapse = "")
}

#Clean printing of data frames with numeric and NA elements.
round_df_char <- function(df, digits, pad = "0", na_vals = "") {
  if (NROW(df) == 0L || NCOL(df) == 0L) {
    return(df)
  }

  if (!is.data.frame(df)) {
    df <- as.data.frame.matrix(df, stringsAsFactors = FALSE)
  }

  rn <- rownames(df)
  cn <- colnames(df)

  infs <- o.negs <- array(FALSE, dim = dim(df))
  nas <- is.na(df)
  nums <- vapply(df, is.numeric, logical(1))

  for (i in which(nums)) {
    infs[, i] <- !nas[, i] & !is.finite(df[[i]])
  }

  for (i in which(!nums)) {
    if (can_str2num(df[[i]])) {
      df[[i]] <- str2num(df[[i]])
      nums[i] <- TRUE
    }
  }

  o.negs[, nums] <- !nas[, nums] & df[nums] < 0 & round(df[nums], digits) == 0
  df[nums] <- round(df[nums], digits = digits)

  pad0 <- identical(as.character(pad), "0")

  for (i in which(nums)) {
    df[[i]] <- format(df[[i]], scientific = FALSE, justify = "none", trim = TRUE,
                      drop0trailing = !pad0)

    if (!pad0 && any(grepl(".", df[[i]], fixed = TRUE))) {
      s <- strsplit(df[[i]], ".", fixed = TRUE)
      lengths <- lengths(s)
      digits.r.of.. <- rep.int(0, NROW(df))
      digits.r.of..[lengths > 1] <- nchar(vapply(s[lengths > 1], `[[`, character(1L), 2))

      dots <- rep.int("", length(s))
      dots[lengths <= 1] <- if (as.character(pad) != "") "." else pad

      pads <- vapply(max(digits.r.of..) - digits.r.of..,
                     function(n) paste(rep.int(pad, n), collapse = ""),
                     character(1L))

      df[[i]] <- paste0(df[[i]], dots, pads)
    }
  }

  df[o.negs] <- paste0("-", df[o.negs])

  # Insert NA placeholders
  df[nas] <- na_vals
  df[infs] <- "N/A"

  if (length(rn) > 0) rownames(df) <- rn
  if (length(cn) > 0) names(df) <- cn

  df
}

#Generalized inverse; port of MASS::ginv()
generalized_inverse <- function(sigma, tol = 1e-8) {
  sigmasvd <- svd(sigma)

  pos <- sigmasvd$d > max(tol * sigmasvd$d[1L], 0)

  sigmasvd$v[, pos, drop = FALSE] %*% (sigmasvd$d[pos]^-1 * t(sigmasvd$u[, pos, drop = FALSE]))
}

#(Weighted) variance that uses special formula for binary variables
wvar <- function(x, bin.var = NULL, w = NULL) {
  if (is_null(w)) w <- rep.int(1, length(x))
  if (is_null(bin.var)) bin.var <- all(x == 0 | x == 1)

  w <- w / sum(w) #weights normalized to sum to 1
  mx <- sum(w * x) #weighted mean

  if (bin.var) {
    return(mx * (1 - mx))
  }

  #Reliability weights variance; same as cov.wt()
  sum(w * (x - mx)^2) / (1 - sum(w^2))
}

#Weighted mean faster than weighted.mean()
wm <- function(x, w = NULL, na.rm = TRUE) {
  if (is_null(w)) {
    if (anyNA(x)) {
      if (!na.rm) return(NA_real_)
      nas <- which(is.na(x))
      x <- x[-nas]
    }
    return(sum(x) / length(x))
  }

  if (anyNA(x) || anyNA(w)) {
    if (!na.rm) return(NA_real_)
    nas <- which(is.na(x) | is.na(w))
    x <- x[-nas]
    w <- w[-nas]
  }

  sum(x * w) / sum(w)
}

#Faster diff()
diff1 <- function(x) {
  x[-1L] - x[-length(x)]
}

#cumsum() for probabilities to ensure they are between 0 and 1
.cumsum_prob <- function(x) {
  s <- cumsum(x)
  s / s[length(s)]
}

#Make vector sum to 1, optionally by group
.make_sum_to_1 <- function(x, by = NULL) {
  if (is_null(by)) {
    return(x / sum(x))
  }

  for (i in unique(by)) {
    in_i <- which(by == i)
    x[in_i] <- x[in_i] / sum(x[in_i])
  }

  x
}

#Make vector sum to n (average of 1), optionally by group
.make_sum_to_n <- function(x, by = NULL) {
  if (is_null(by)) {
    return(length(x) * x / sum(x))
  }

  for (i in unique(by)) {
    in_i <- which(by == i)
    x[in_i] <- length(in_i) * x[in_i] / sum(x[in_i])
  }

  x
}

#Generate a new list, data.frame, or matrix
make_list <- function(n) {
  if (is_null(n)) {
    vector("list", 0L)
  }
  else if (rlang::is_scalar_integerish(n) && n >= 0) {
    vector("list", as.integer(n))
  }
  else if (is.atomic(n)) {
    setNames(vector("list", length(n)),
             as.character(n))
  }
  else {
    stop("'n' must be an integer(ish) scalar or an atomic variable.")
  }
}
make_df <- function(ncol, nrow = 0L, types = "numeric") {
  if (missing(ncol) || is_null(ncol)) {
    ncol <- 0L
  }

  if (rlang::is_scalar_integerish(ncol)) {
    col_names <- NULL
    ncol <- as.integer(ncol)
  }
  else if (is.atomic(ncol)) {
    col_names <- as.character(ncol)
    ncol <- length(ncol)
  }

  if (is_null(nrow)) {
    nrow <- 0L
  }

  if (rlang::is_scalar_integerish(nrow)) {
    row_names <- NULL
    nrow <- as.integer(nrow)
  }
  else if (is.atomic(nrow)) {
    row_names <- as.character(nrow)
    nrow <- length(nrow)
  }

  df <- as.data.frame.matrix(matrix(NA_real_, nrow = nrow, ncol = ncol))

  names(df) <- col_names
  rownames(df) <- row_names

  if (is_null(types)) {
    return(df)
  }

  if (!length(types) %in% c(1L, ncol)) {
    stop("'types' must be equal to the number of columns.")
  }

  if (!is.character(types) ||
      !all(types %in% c("numeric", "integer", "logical", "character", NA))) {
    stop("'types' must be an acceptable type. For factors, use NA.")
  }

  if (length(types) == 1L) {
    types <- rep.int(types, ncol)
  }

  for (i in which(!is.na(types))) {
    df[[i]][] <- switch(types[i],
                        integer = NA_integer_,
                        logical = NA,
                        character = NA_character_,
                        NA_real_)
  }

  df
}
make_matrix <- function(ncol, nrow = 0L, type = "numeric") {
  if (missing(ncol) || is_null(ncol)) {
    ncol <- 0L
  }

  if (rlang::is_scalar_integerish(ncol)) {
    col_names <- NULL
    ncol <- as.integer(ncol)
  }
  else if (is.atomic(ncol)) {
    col_names <- as.character(ncol)
    ncol <- length(ncol)
  }

  if (is_null(nrow)) {
    nrow <- 0L
  }

  if (rlang::is_scalar_integerish(nrow)) {
    row_names <- NULL
    nrow <- as.integer(nrow)
  }
  else if (is.atomic(nrow)) {
    row_names <- as.character(nrow)
    nrow <- length(nrow)
  }

  if (is_null(type)) {
    type <- NA
  }

  if (length(type) != 1L) {
    stop("'type' must have length 1.")
  }

  if (!is.character(type) ||
      (!is.na(type) && !(type %in% c("numeric", "integer", "logical", "character")))) {
    stop("'types' must be an acceptable type. For factors, use NA.")
  }

  fill_val <- switch(type,
                     integer = NA_integer_,
                     logical = NA,
                     character = NA_character_,
                     NA_real_)

  matrix(fill_val, nrow = nrow, ncol = ncol,
         dimnames = list(row_names, col_names))
}

#Extract variables from ..., similar to ...elt() or get0(), by name without evaluating list(...)
...get <- function(x, ifnotfound = NULL) {
  expr <- quote({
    .m1 <- match(.x, ...names())
    if (anyNA(.m1)) {
      .ifnotfound
    }
    else {
      .m2 <- ...elt(.m1[1L])
      if (is_not_null(.m2)) .m2
      else .ifnotfound
    }
  })

  eval(expr,
       pairlist(.x = x[1L], .ifnotfound = ifnotfound),
       parent.frame(1L))
}

#Extract multiple variables from ..., similar to mget(), by name without evaluating list(...)
...mget <- function(x) {
  found <- match(x, eval(quote(...names()), parent.frame(1L)))

  not_found <- is.na(found)

  if (all(not_found)) {
    return(list())
  }

  setNames(lapply(found[!not_found], function(z) {
    eval(quote(...elt(.z)),
         pairlist(.z = z),
         parent.frame(3L))
  }), x[!not_found])
}

#Helper function to fill named vectors with x and given names of y
rep_with <- function(x, y) {
  setNames(rep.int(x, length(y)), names(y))
}

#cat() if verbose = TRUE (default sep = "", line wrapping)
.cat_verbose <- function(..., verbose = TRUE, sep = "") {
  if (!verbose) {
    return(invisible(NULL))
  }

  m <- do.call(function(...) paste(..., sep = sep), list(...))

  if (endsWith(m, "\n")) {
    m <- paste0(paste(strwrap(m), collapse = "\n"), "\n")
  }
  else {
    m <- paste(strwrap(m), collapse = "\n")
  }

  cat(paste(m, collapse = "\n"))
}
