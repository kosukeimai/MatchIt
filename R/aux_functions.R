#Auxiliary functions; some from WeightIt

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
    stop(sprintf("Not enough units to fit treated and control %s in each subclass.",
                 ngettext(min.n, "unit", "units")), call. = FALSE)
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
create_info <- function(method, fn1, link, discard, replace, ratio, mahalanobis, transform, subclass, antiexact, distance_is_matrix) {
  info <- list(method = method,
               distance = if (is.null(fn1)) NULL else sub("distance2", "", fn1, fixed = TRUE),
               link = if (is.null(link)) NULL else link,
               discard = discard,
               replace = if (!is.null(method) && method %in% c("nearest", "genetic")) replace else NULL,
               ratio = if (!is.null(method) && method %in% c("nearest", "optimal", "genetic")) ratio else NULL,
               max.controls = if (!is.null(method) && method %in% c("nearest", "optimal")) attr(ratio, "max.controls") else NULL,
               mahalanobis = mahalanobis,
               transform = transform,
               subclass = if (!is.null(method) && method == "subclass") length(unique(subclass[!is.na(subclass)])) else NULL,
               antiexact = antiexact,
               distance_is_matrix = distance_is_matrix)
  info
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
           "quick" = "generalized full matching",
           "genetic" = "genetic matching",
           "subclass" = paste0("subclassification (", info$subclass, " subclasses)"),
           "cardinality" = "cardinality matching",
           if (is.null(attr(info$method, "method"))) "an unspecified matching method"
           else attr(info$method, "method"))
  out.list[["replace"]] <- if (!is.null(info$replace) && info$method %in% c("nearest", "genetic")) {
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
    else dist <- sprintf("GLM with a %s link", link)
  }
  else if (distance == "gam") {
    dist <- sprintf("GAM with a %s link", link)
  }
  else if (distance == "gbm") {
    dist <- "GBM"
  }
  else if (distance == "elasticnet") {
    dist <- sprintf("an elastic net with a %s link", link)
  }
  else if (distance == "lasso") {
    if (link == "logit") dist <- "lasso logistic regression"
    else dist <- sprintf("lasso regression with a %s link", link)
  }
  else if (distance == "ridge") {
    dist <- sprintf("ridge regression with a %s link", link)
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
    attr(out, "plural") <- FALSE
  }
  else {
    word.list <- word.list[!word.list %in% c(NA_character_, "")]
    L <- length(word.list)
    if (L == 0) {
      out <- ""
      attr(out, "plural") <- FALSE
    }
    else if (L == 1) {
      out <- word.list
      if (is.are) out <- paste(out, "is")
      attr(out, "plural") <- FALSE
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
      attr(out, "plural") <- TRUE
    }

  }
  return(out)
}

#Add quotation marks around a string.
add_quotes <- function(x, quotes = 2L) {
  if (!isFALSE(quotes)) {
    if (isTRUE(quotes) || as.integer(quotes) == 2L) x <- sprintf('"%s"', x)
    else if (as.integer(quotes) == 1L) x <- sprintf("'%s'", x)
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
  arg.name <- deparse1(substitute(arg), width.cutoff = 500L)

  if (missing(choices)) {
    formal.args <- formals(sys.function(sysP <- sys.parent()))
    choices <- eval(formal.args[[as.character(substitute(arg))]],
                    envir = sys.frame(sysP))
  }

  if (is.null(arg))
    return(choices[1L])
  else if (!is.character(arg))
    stop(sprintf("The argument to '%s' must be NULL or a character vector", arg.name), call. = FALSE)
  if (!several.ok) {
    if (identical(arg, choices))
      return(arg[1L])
    if (length(arg) > 1L)
      stop(sprintf("The argument to '%s' must be of length 1", arg.name), call. = FALSE)
  }
  else if (length(arg) == 0)
    stop(sprintf("The argument to '%s' must be of length >= 1", arg.name), call. = FALSE)

  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
  if (all(i == 0L))
    stop(sprintf("The argument to '%s' should be %s %s.",
                arg.name, ngettext(length(choices), "", if (several.ok) "at least one of " else "one of "),
                word_list(choices, and.or = "or", quotes = 2)),
         call. = FALSE)
  i <- i[i > 0L]
  if (!several.ok && length(i) > 1)
    stop("There is more than one match in 'match_arg'")
  choices[i]
}

#Turn a vector into a 0/1 vector. 'zero' and 'one' can be supplied to make it clear which is
#which; otherwise, a guess is used. From WeightIt.
binarize <- function(variable, zero = NULL, one = NULL) {
  var.name <- deparse1(substitute(variable))
  if (length(unique(variable)) > 2) {
    stop(sprintf("Cannot binarize %s: more than two levels.", var.name), call. = FALSE)
  }
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

#Make interaction vector out of matrix of covs; similar to interaction()
exactify <- function(X, nam = NULL, sep = "|", include_vars = FALSE) {
  if (is.null(nam)) nam <- rownames(X)
  if (is.matrix(X)) X <- setNames(lapply(seq_len(ncol(X)), function(i) X[,i]), colnames(X))
  if (!is.list(X)) stop("X must be a matrix, data frame, or list.")

  if (include_vars) {
    for (i in seq_along(X)) {
      if (is.character(X[[i]]) || is.factor(X[[i]])) {
        X[[i]] <- sprintf('%s = "%s"', names(X)[i], X[[i]])
      }
      else {
        X[[i]] <- sprintf('%s = %s', names(X)[i], X[[i]])
      }
    }
  }
  else {
    for (i in seq_along(X)) {
      if (is.factor(X[[i]])) {
        X[[i]] <- format(levels(X[[i]]), justify = "right")[X[[i]]]
      }
      else {
        X[[i]] <- format(X[[i]], justify = "right")
      }
    }
  }

  out <- do.call("paste", c(X, sep = sep))
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
  is.na(x_num)[nas] <- TRUE
  return(x_num)
}

#Capitalize first letter of string
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

#Capitalize first letter of each word
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste0(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           collapse = " ")
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

#Clean printing of data frames with numeric and NA elements.
round_df_char <- function(df, digits, pad = "0", na_vals = "") {
  #Digits is passed to round(). pad is used to replace trailing zeros so decimal
  #lines up. Should be "0" or " "; "" (the empty string) un-aligns decimals.
  #na_vals is what NA should print as.

  if (!is.data.frame(df)) df <- as.data.frame.matrix(df, stringsAsFactors = FALSE)
  if (NROW(df) == 0 || NCOL(df) == 0) return(df)
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

  df <- as.matrix(df)

  return(df)
}

#Generalized inverse; port of MASS::ginv()
generalized_inverse <- function(sigma) {
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

  mf <- droplevels(mf)

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
  charmm <- array(NA_character_, dim = dim(nummm), dimnames = dimnames(nummm))
  charmm[] <- names(treat)[nummm]
  charmm
}

charmm2nummm <- function(charmm, treat) {
  nummm <- array(NA_integer_, dim = dim(charmm))
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

#Pooled within-group (weighted) covariance by group-mean centering covariates. Used
#in Mahalanobis distance
pooled_cov <- function(X, t, w = NULL) {
  unique_t <- unique(t)
  if (is.null(dim(X))) X <- matrix(X, nrow = length(X))
  if (is.null(w)) {
    n <- nrow(X)
    for (i in unique_t) {
      in_t <- which(t == i)
      for (j in seq_len(ncol(X))) {
        X[in_t, j] <- X[in_t, j] - mean(X[in_t, j])
      }
    }
    pooled_cov <- cov(X)*(n-1)/(n-length(unique_t))
  }
  else {
    for (i in unique_t) {
      in_t <- which(t == i)
      for (j in seq_len(ncol(X))) {
        X[in_t, j] <- X[in_t, j] - wm(X[in_t, j], w[in_t])
      }
    }
    pooled_cov <- cov.wt(X, w)$cov
  }
  return(pooled_cov)
}

pooled_sd <- function(X, t, w = NULL, bin.var = NULL, contribution = "proportional") {
  contribution <- match_arg(contribution, c("proportional", "equal"))
  unique_t <- unique(t)
  if (is.null(dim(X))) X <- matrix(X, nrow = length(X))
  n <- nrow(X)
  if (is.null(bin.var)) bin.var <- apply(X, 2, function(x) all(x == 0 | x == 1))

  if (contribution == "equal") {
    vars <- matrix(0, nrow = length(unique_t), ncol = ncol(X))
    for (i in seq_along(unique_t)) {
      in_t <- which(t == unique_t[i])
      vars[i,] <- vapply(seq_len(ncol(X)), function(j) {
        x <- X[,j]
        b <- bin.var[j]
        wvar(x[in_t], w = w[in_t], bin.var = b)
      }, numeric(1L))
    }
    pooled_var <- colMeans(vars)
  }
  else {
    pooled_var <- vapply(seq_len(ncol(X)), function(j) {
      x <- X[,j]
      b <- bin.var[j]

      if (b) {
        if (is.null(w)) {
          v <- vapply(unique_t, function(i) {
            sxi <- sum(x[t == i])
            ni <- sum(t == i)
            sxi * (1 - sxi/ni) / n
          }, numeric(1L))
          return(sum(v))
        }
        else {
          v <- vapply(unique_t, function(i) {
            sxi <- sum(x[t == i] * w[t == i])
            ni <- sum(w[t==i])
            sxi * (1 - sxi/ni) / sum(w)
          }, numeric(1L))
          return(sum(v))
        }
      }
      else {
        if (is.null(w)) {
          for (i in unique_t) {
            x[t==i] <- x[t==i] - wm(x[t==i])
          }
          return(sum(x^2)/(n - length(unique_t)))
        }
        else {
          for (i in unique_t) {
            x[t==i] <- x[t==i] - wm(x[t==i], w[t==i])
          }
          w_ <- w/sum(w)
          return(sum(w_ * x^2)/(1 - sum(w_^2)))
        }
      }
    }, numeric(1L))
  }
  return(sqrt(pooled_var))
}

#Effective sample size
ESS <- function(w) {
  sum(w)^2/sum(w^2)
}

#Compute sample sizes
nn <- function(treat, weights, discarded = NULL, s.weights = NULL) {

  if (is.null(discarded)) discarded <- rep(FALSE, length(treat))
  if (is.null(s.weights)) s.weights <- rep(1, length(treat))
  weights <- weights * s.weights
  n <- matrix(0, ncol=2, nrow=6, dimnames = list(c("All (ESS)", "All", "Matched (ESS)","Matched", "Unmatched","Discarded"),
                                                 c("Control", "Treated")))

  #                      Control                                    Treated
  n["All (ESS)",] <-     c(ESS(s.weights[treat==0]),                ESS(s.weights[treat==1]))
  n["All",] <-           c(sum(treat==0),                           sum(treat==1))
  n["Matched (ESS)",] <- c(ESS(weights[treat==0]),                  ESS(weights[treat==1]))
  n["Matched",] <-       c(sum(treat==0 & weights > 0),             sum(treat==1 & weights > 0))
  n["Unmatched",] <-     c(sum(treat==0 & weights==0 & !discarded), sum(treat==1 & weights==0 & !discarded))
  n["Discarded",] <-     c(sum(treat==0 & discarded),               sum(treat==1 & discarded))

  return(n)
}

#Compute subclass sample sizes
qn <- function(treat, subclass, discarded = NULL) {

  treat <- factor(treat, levels = 0:1, labels = c("Control", "Treated"))
  if (is.null(discarded)) discarded <- rep(FALSE, length(treat))
  qn <- table(treat[!discarded], subclass[!discarded])

  if (any(is.na(subclass) & !discarded)) {
    qn <- cbind(qn, table(treat[is.na(subclass) & !discarded]))
    colnames(qn)[ncol(qn)] <- "Unmatched"
  }
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

#Faster diff()
diff1 <- function(x) {
  x[-1] - x[-length(x)]
}

#Check if is a whole number
is_whole_number <- function(x) {
  abs(x - round(x)) < .Machine$double.eps
}
