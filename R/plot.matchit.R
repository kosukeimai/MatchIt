#' Generate Balance Plots after Matching and Subclassification
#'
#' Generates plots displaying distributional balance and overlap on covariates
#' and propensity scores before and after matching and subclassification. For
#' displaying balance solely on covariate standardized mean differences, see
#' [plot.summary.matchit()]. The plots here can be used to assess to what
#' degree covariate and propensity score distributions are balanced and how
#' weighting and discarding affect the distribution of propensity scores.
#'
#' @aliases plot.matchit plot.matchit.subclass
#'
#' @param x a `matchit` object; the output of a call to [matchit()].
#' @param type the type of plot to display. Options include `"qq"`,
#' `"ecdf"`, `"density"`, `"jitter"`, and `"histogram"`.
#' See Details. Default is `"qq"`. Abbreviations allowed.
#' @param interactive `logical`; whether the graphs should be displayed in
#' an interactive way. Only applies for `type = "qq"`, `"ecdf"`,
#' `"density"`, and `"jitter"`. See Details.
#' @param which.xs with `type = "qq"`, `"ecdf"`, or `"density"`,
#' for which covariate(s) plots should be displayed. Factor variables should be
#' named by the original variable name rather than the names of individual
#' dummy variables created after expansion with `model.matrix`. Can be supplied as a character vector or a one-sided formula.
#' @param data an optional data frame containing variables named in `which.xs` but not present in the `matchit` object.
#' @param subclass with subclassification and `type = "qq"`,
#' `"ecdf"`, or `"density"`, whether to display balance for
#' individual subclasses, and, if so, for which ones. Can be `TRUE`
#' (display plots for all subclasses), `FALSE` (display plots only in
#' aggregate), or the indices (e.g., `1:6`) of the specific subclasses for
#' which to display balance. When unspecified, if `interactive = TRUE`,
#' you will be asked for which subclasses plots are desired, and otherwise,
#' plots will be displayed only in aggregate.
#' @param \dots arguments passed to [plot()] to control the appearance of the
#' plot. Not all options are accepted.
#'
#' @details
#' `plot.matchit()` makes one of five different plots depending on the
#' argument supplied to `type`. The first three, `"qq"`,
#' `"ecdf"`, and `"density"`, assess balance on the covariates. When
#' `interactive = TRUE`, plots for three variables will be displayed at a
#' time, and the prompt in the console allows you to move on to the next set of
#' variables. When `interactive = FALSE`, multiple pages are plotted at
#' the same time, but only the last few variables will be visible in the
#' displayed plot. To see only a few specific variables at a time, use the
#' `which.xs` argument to display plots for just those variables. If fewer
#' than three variables are available (after expanding factors into their
#' dummies), `interactive` is ignored.
#'
#' With `type = "qq"`, empirical quantile-quantile (eQQ) plots are created
#' for each covariate before and after matching. The plots involve
#' interpolating points in the smaller group based on the weighted quantiles of
#' the other group. When points are approximately on the 45-degree line, the
#' distributions in the treatment and control groups are approximately equal.
#' Major deviations indicate departures from distributional balance. With
#' variable with fewer than 5 unique values, points are jittered to more easily
#' visualize counts.
#'
#' With `type = "ecdf"`, empirical cumulative density function (eCDF)
#' plots are created for each covariate before and after matching. Two eCDF
#' lines are produced in each plot: a gray one for control units and a black
#' one for treated units. Each point on the lines corresponds to the proportion
#' of units (or proportionate share of weights) less than or equal to the
#' corresponding covariate value (on the x-axis). Deviations between the lines
#' on the same plot indicates distributional imbalance between the treatment
#' groups for the covariate. The eCDF and eQQ statistics in [summary.matchit()]
#' correspond to these plots: the eCDF max (also known as the
#' Kolmogorov-Smirnov statistic) and mean are the largest and average vertical
#' distance between the lines, and the eQQ max and mean are the largest and
#' average horizontal distance between the lines.
#'
#' With `type = "density"`, density plots are created for each covariate
#' before and after matching. Two densities are produced in each plot: a gray
#' one for control units and a black one for treated units. The x-axis
#' corresponds to the value of the covariate and the y-axis corresponds to the
#' density or probability of that covariate value in the corresponding group.
#' For binary covariates, bar plots are produced, having the same
#' interpretation. Deviations between the black and gray lines represent
#' imbalances in the covariate distribution; when the lines coincide (i.e.,
#' when only the black line is visible), the distributions are identical.
#'
#' The last two plots, `"jitter"` and `"histogram"`, visualize the
#' distance (i.e., propensity score) distributions. These plots are more for
#' heuristic purposes since the purpose of matching is to achieve balance on
#' the covariates themselves, not the propensity score.
#'
#' With `type = "jitter"`, a jitter plot is displayed for distance values
#' before and after matching. This method requires a distance variable (e.g., a
#' propensity score) to have been estimated or supplied in the call to
#' `matchit()`. The plot displays individuals values for matched and
#' unmatched treatment and control units arranged horizontally by their
#' propensity scores. Points are jitter so counts are easier to see. The size
#' of the points increases when they receive higher weights. When
#' `interactive = TRUE`, you can click on points in the graph to identify
#' their rownames and indices to further probe extreme values, for example.
#' With subclassification, vertical lines representing the subclass boundaries
#' are overlay on the plots.
#'
#' With `type = "histogram"`, a histogram of distance values is displayed
#' for the treatment and control groups before and after matching. This method
#' requires a distance variable (e.g., a propensity score) to have been
#' estimated or supplied in the call to `matchit()`. With
#' subclassification, vertical lines representing the subclass boundaries are
#' overlay on the plots.
#'
#' With all methods, sampling weights are incorporated into the weights if
#' present.
#'
#' @note Sometimes, bugs in the plotting functions can cause strange layout or
#' size issues. Running [frame()] or [dev.off()] can be used to reset the
#' plotting pane (note the latter will delete any plots in the plot history).
#'
#' @seealso [summary.matchit()] for numerical summaries of balance, including
#' those that rely on the eQQ and eCDF plots.
#'
#' [plot.summary.matchit()] for plotting standardized mean differences in a
#' Love plot.
#'
#' \pkgfun{cobalt}{bal.plot} for displaying distributional balance in several other
#' ways that are more easily customizable and produce *ggplot2* objects.
#' *cobalt* functions natively support `matchit` objects.
#'
#' @examples
#' data("lalonde")
#'
#' m.out <- matchit(treat ~ age + educ + married +
#'                    race + re74, data = lalonde,
#'                  method = "nearest")
#' plot(m.out, type = "qq", interactive = FALSE,
#'      which.xs = ~age + educ + re74)
#' plot(m.out, type = "histogram")
#'
#' s.out <- matchit(treat ~ age + educ + married +
#'                    race + nodegree + re74 + re75,
#'                  data = lalonde, method = "subclass")
#' plot(s.out, type = "density", interactive = FALSE,
#'      which.xs = ~age + educ + re74,
#'      subclass = 3)
#' plot(s.out, type = "jitter", interactive = FALSE)
#'

#' @exportS3Method plot matchit
plot.matchit <- function(x, type = "qq", interactive = TRUE, which.xs = NULL, data = NULL, ...) {

  chk::chk_string(type)
  type <- tolower(type)
  type <- match_arg(type, c("qq", "ecdf", "density", "jitter", "histogram"))

  if (type %in% c("qq", "ecdf", "density")) {
    matchit.covplot(x, type = type, interactive = interactive,
                    which.xs = which.xs, data = data, ...)
  }
  else if (type == "jitter") {
    if (is.null(x$distance)) {
      .err("`type = \"jitter\"` cannot be used if a distance measure is not estimated or supplied. No plots generated")
    }
    jitter.pscore(x, interactive = interactive,...)
  }
  else if (type == "histogram") {
    if (is.null(x$distance)) {
      .err("`type = \"hist\"` cannot be used if a distance measure is not estimated or supplied. No plots generated")
    }
    hist.pscore(x,...)
  }
  invisible(x)
}

#' @exportS3Method plot matchit.subclass
#' @rdname plot.matchit
plot.matchit.subclass <- function(x, type = "qq", interactive = TRUE, which.xs = NULL, subclass, ...) {
  choice.menu <- function(choices, question) {
    k <- length(choices)-1
    Choices <- data.frame(choices)
    row.names(Choices) <- 0:k
    names(Choices) <- "Choices"
    print.data.frame(Choices, right=FALSE)
    ans <- readline(question)
    while (!ans %in% 0:k) {
      message("Not valid -- please pick one of the choices")
      print.data.frame(Choices, right=FALSE)
      ans <- readline(question)
    }

    ans
  }

  type <- tolower(type)
  type <- match_arg(type, c("qq", "ecdf", "density", "jitter", "histogram"))

  if (type %in% c("qq", "ecdf", "density")) {
    #If subclass = T, index, or range, display all or range of subclasses, using interactive to advance
    #If subclass = F, display aggregate across subclass, using interactive to advance
    #If subclass = NULL, if interactive, use to choose subclass, else display aggregate across subclass

    subclasses <- levels(x$subclass)
    miss.sub <- missing(subclass) || is.null(subclass)
    if (miss.sub || isFALSE(subclass)) which.subclass <- NULL
    else if (isTRUE(subclass)) which.subclass <- subclasses
    else if (!is.atomic(subclass) || !all(subclass %in% seq_along(subclasses))) {
      .err("`subclass` should be `TRUE`, `FALSE`, or a vector of subclass indices for which subclass balance is to be displayed")
    }
    else which.subclass <- subclasses[subclass]

    if (!is.null(which.subclass)) {
      matchit.covplot.subclass(x, type = type, which.subclass = which.subclass,
                               interactive = interactive, which.xs = which.xs, ...)
    }
    else if (interactive && miss.sub) {
      subclasses <- levels(x$subclass)
      choices <- c("No (Exit)", paste0("Yes: Subclass ", subclasses), "Yes: In aggregate")
      plot.name <- switch(type, "qq" = "quantile-quantile", "ecdf" = "empirical CDF", "density" = "density")
      question <- sprintf("Would you like to see %s plots of any subclasses? ", plot.name)
      ans <- -1
      while (ans != 0) {
        ans <- as.numeric(choice.menu(choices, question))
        if (ans %in% seq_along(subclasses) && any(x$subclass == subclasses[ans])) {
          matchit.covplot.subclass(x, type = type, which.subclass = subclasses[ans],
                                   interactive = interactive, which.xs = which.xs, ...)
        }
        else if (ans != 0) {
          matchit.covplot(x, type = type, interactive = interactive, which.xs = which.xs, ...)
        }
      }
    }
    else {
      matchit.covplot(x, type = type, interactive = interactive, which.xs = which.xs, ...)
    }
  }
  else if (type=="jitter") {
    if (is.null(x$distance)) {
      .err("`type = \"jitter\"` cannot be used when no distance variable was estimated or supplied")
    }
    jitter.pscore(x, interactive = interactive, ...)
  }
  else if (type == "histogram") {
    if (is.null(x$distance)) {
      .err("`type = \"histogram\"` cannot be used when no distance variable was estimated or supplied")
    }
    hist.pscore(x,...)
  }
  invisible(x)
}

## plot helper functions
matchit.covplot <- function(object, type = "qq", interactive = TRUE, which.xs = NULL, data = NULL, ...) {

  if (is.null(which.xs)) {
    if (length(object$X) == 0) {
      .wrn("No covariates to plot")
      return(invisible(NULL))
    }
    X <- object$X

    if (!is.null(object$exact)) {
      Xexact <- model.frame(object$exact, data = object$X)
      X <- cbind(X, Xexact[setdiff(names(Xexact), names(X))])
    }

    if (!is.null(object$mahvars)) {
      Xmahvars <- model.frame(object$mahvars, data = object$X)
      X <- cbind(X, Xmahvars[setdiff(names(Xmahvars), names(X))])
    }
  }
  else {
    if (!is.null(data)) {
      if (!is.data.frame(data) || nrow(data) != length(object$treat)) {
        .err("`data` must be a data frame with as many rows as there are units in the supplied `matchit` object")
      }
      data <- cbind(data, object$X[setdiff(names(object$X), names(data))])
    }
    else {
      data <- object$X
    }

    if (!is.null(object$exact)) {
      Xexact <- model.frame(object$exact, data = object$X)
      data <- cbind(data, Xexact[setdiff(names(Xexact), names(data))])
    }

    if (!is.null(object$mahvars)) {
      Xmahvars <- model.frame(object$mahvars, data = object$X)
      data <- cbind(data, Xmahvars[setdiff(names(Xmahvars), names(data))])
    }

    if (is.character(which.xs)) {
      if (!all(which.xs %in% names(data))) {
        .err("All variables in `which.xs` must be in the supplied `matchit` object or in `data`")
      }
      X <- data[which.xs]
    }
    else if (rlang::is_formula(which.xs)) {
      which.xs <- update(which.xs, NULL ~ .)
      X <- model.frame(which.xs, data, na.action = "na.pass")
    }
    else {
      .err("`which.xs` must be supplied as a character vector of names or a one-sided formula")
    }

    # if (anyNA(X)) {
    #   stop("Missing values are not allowed in the covariates named in `which.xs`.",
    #        call. = FALSE)
    # }

    k <- ncol(X)
    for (i in seq_len(k)) {
      if (anyNA(X[[i]]) || (is.numeric(X[[i]]) && any(!is.finite(X[[i]])))) {
        covariates.with.missingness <- names(X)[i:k][vapply(i:k, function(j) anyNA(X[[j]]) ||
                                                                          (is.numeric(X[[j]]) &&
                                                                             any(!is.finite(X[[j]]))),
                                                                        logical(1L))]
        .err(paste0("Missing and non-finite values are not allowed in the covariates named in `which.xs`. Variables with missingness or non-finite values:\n\t",
                    paste(covariates.with.missingness, collapse = ", ")), tidy = FALSE)
      }
      if (is.character(X[[i]])) X[[i]] <- factor(X[[i]])
    }
  }

  # chars.in.X <- vapply(X, is.character, logical(1L))
  # X[chars.in.X] <- lapply(X[chars.in.X], factor)

  X <- droplevels(X)

  t <- object$treat

  sw <- if (is.null(object$s.weights)) rep(1, length(t)) else object$s.weights
  w <- object$weights * sw
  if (is.null(w)) w <- rep(1, length(t))

  w <- .make_sum_to_1(w, by = t)
  sw <- .make_sum_to_1(sw, by = t)

  if (type == "density") {
    varnames <- names(X)
  }
  else {
    X <- get.covs.matrix(data = X)
    varnames <- colnames(X)
  }

  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  oma <- c(2.25, 0, 3.75, 1.5)
  if (type == "qq") {
    opar <- par(mfrow = c(3, 3), mar = rep.int(1/2, 4), oma = oma)
  }
  else if (type %in% c("ecdf", "density")) {
    opar <- par(mfrow = c(3, 3), mar = c(1.5,.5,1.5,.5), oma = oma)
  }

  for (i in seq_along(varnames)){
    x <- if (type == "density") X[[i]] else X[,i]

    plot.new()

    if (((i-1)%%3)==0) {

      if (type == "qq") {
        htext <- "eQQ Plots"
        mtext(htext, 3, 2, TRUE, 0.5, cex=1.1, font = 2)
        mtext("All", 3, .25, TRUE, 0.5, cex=1, font = 1)
        mtext("Matched", 3, .25, TRUE, 0.83, cex=1, font = 1)
        mtext("Control Units", 1, 0, TRUE, 2/3, cex=1, font = 1)
        mtext("Treated Units", 4, 0, TRUE, 0.5, cex=1, font = 1)
      }
      else if (type == "ecdf") {
        htext <- "eCDF Plots"
        mtext(htext, 3, 2, TRUE, 0.5, cex=1.1, font = 2)
        mtext("All", 3, .25, TRUE, 0.5, cex=1, font = 1)
        mtext("Matched", 3, .25, TRUE, 0.83, cex=1, font = 1)
      }
      else if (type == "density") {
        htext <- "Density Plots"
        mtext(htext, 3, 2, TRUE, 0.5, cex=1.1, font = 2)
        mtext("All", 3, .25, TRUE, 0.5, cex=1, font = 1)
        mtext("Matched", 3, .25, TRUE, 0.83, cex=1, font = 1)
      }

    }

    par(usr = c(0, 1, 0, 1))
    l.wid <- strwidth(varnames, "user")
    cex.labels <- max(0.75, min(1.45, 0.85/max(l.wid)))
    text(0.5, 0.5, varnames[i], cex = cex.labels)

    if (type == "qq") {
      qqplot_match(x = x, t = t, w = w, sw = sw, ...)
    }
    else if (type == "ecdf") {
      ecdfplot_match(x = x, t = t, w = w, sw = sw, ...)
    }
    else if (type == "density") {
      densityplot_match(x = x, t = t, w = w, sw = sw, ...)
    }

    devAskNewPage(ask = interactive)
  }
  devAskNewPage(ask = FALSE)

  invisible(NULL)
}

matchit.covplot.subclass <- function(object, type = "qq", which.subclass = NULL,
                                     interactive = TRUE, which.xs = NULL, data = NULL, ...) {

  if (is.null(which.xs)) {
    if (length(object$X) == 0) {
      .wrn("No covariates to plot")
      return(invisible(NULL))
    }
    X <- object$X

    if (!is.null(object$exact)) {
      Xexact <- model.frame(object$exact, data = object$X)
      X <- cbind(X, Xexact[setdiff(names(Xexact), names(X))])
    }

    if (!is.null(object$mahvars)) {
      Xmahvars <- model.frame(object$mahvars, data = object$X)
      X <- cbind(X, Xmahvars[setdiff(names(Xmahvars), names(X))])
    }
  }
  else {
    if (!is.null(data)) {
      if (!is.data.frame(data) || nrow(data) != length(object$treat)) {
        .err("`data` must be a data frame with as many rows as there are units in the supplied `matchit` object")
      }
      data <- cbind(data, object$X[setdiff(names(object$X), names(data))])
    }
    else {
      data <- object$X
    }

    if (!is.null(object$exact)) {
      Xexact <- model.frame(object$exact, data = object$X)
      data <- cbind(data, Xexact[setdiff(names(Xexact), names(data))])
    }

    if (!is.null(object$mahvars)) {
      Xmahvars <- model.frame(object$mahvars, data = object$X)
      data <- cbind(data, Xmahvars[setdiff(names(Xmahvars), names(data))])
    }

    if (is.character(which.xs)) {
      if (!all(which.xs %in% names(data))) {
        .err("All variables in `which.xs` must be in the supplied `matchit` object or in `data`")
      }
      X <- data[which.xs]
    }
    else if (inherits(which.xs, "formula")) {
      which.xs <- update(which.xs, NULL ~ .)
      X <- model.frame(which.xs, data, na.action = "na.pass")

      if (anyNA(X)) {
        .err("Missing values are not allowed in the covariates named in `which.xs`")
      }
    }
    else {
      .err("`which.xs` must be supplied as a character vector of names or a one-sided formula")
    }
  }

  chars.in.X <- vapply(X, is.character, logical(1L))
  X[chars.in.X] <- lapply(X[chars.in.X], factor)

  X <- droplevels(X)

  t <- object$treat

  if (!is.atomic(which.subclass)) {
    .err("The argument to `subclass` must be NULL or the indices of the subclasses for which to display covariate distributions")
  }
  if (!all(which.subclass %in% object$subclass[!is.na(object$subclass)])) {
    .err("The argument supplied to `subclass` is not the index of any subclass in the matchit object")
  }

  if (type == "density") {
    varnames <- names(X)
  }
  else {
    X <- get.covs.matrix(data = X)
    varnames <- colnames(X)
  }

  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  oma <- c(2.25, 0, 3.75, 1.5)

  for (s in which.subclass) {
    if (type == "qq") {
      opar <- par(mfrow = c(3, 3), mar = rep.int(1/2, 4), oma = oma)
    }
    else if (type %in% c("ecdf", "density")) {
      opar <- par(mfrow = c(3, 3), mar = c(1.5,.5,1.5,.5), oma = oma)
    }

    sw <- if (is.null(object$s.weights)) rep(1, length(t)) else object$s.weights
    w <- sw*(!is.na(object$subclass) & object$subclass == s)

    w <- .make_sum_to_1(w, by = t)
    sw <- .make_sum_to_1(sw, by = t)

    for (i in seq_along(varnames)){

      x <- if (type == "density") X[[i]] else X[,i]

      plot.new()

      if (((i-1)%%3)==0) {

        if (type == "qq") {
          htext <- paste0("eQQ Plots (Subclass ", s,")")
          mtext(htext, 3, 2, TRUE, 0.5, cex=1.1, font = 2)
          mtext("All", 3, .25, TRUE, 0.5, cex=1, font = 1)
          mtext("Matched", 3, .25, TRUE, 0.83, cex=1, font = 1)
          mtext("Control Units", 1, 0, TRUE, 2/3, cex=1, font = 1)
          mtext("Treated Units", 4, 0, TRUE, 0.5, cex=1, font = 1)
        }
        else if (type == "ecdf") {
          htext <- paste0("eCDF Plots (Subclass ", s,")")
          mtext(htext, 3, 2, TRUE, 0.5, cex=1.1, font = 2)
          mtext("All", 3, .25, TRUE, 0.5, cex=1, font = 1)
          mtext("Matched", 3, .25, TRUE, 0.83, cex=1, font = 1)
        }
        else if (type == "density") {
          htext <- paste0("Density Plots (Subclass ", s,")")
          mtext(htext, 3, 2, TRUE, 0.5, cex=1.1, font = 2)
          mtext("All", 3, .25, TRUE, 0.5, cex=1, font = 1)
          mtext("Matched", 3, .25, TRUE, 0.83, cex=1, font = 1)
        }

      }

      #Empty plot with variable name
      par(usr = c(0, 1, 0, 1))
      l.wid <- strwidth(varnames, "user")
      cex.labels <- max(0.75, min(1.45, 0.85/max(l.wid)))
      text(0.5, 0.5, varnames[i], cex = cex.labels)

      if (type == "qq") {
        qqplot_match(x = x, t = t, w = w, sw = sw, ...)
      }
      else if (type == "ecdf") {
        ecdfplot_match(x = x, t = t, w = w, sw = sw, ...)
      }
      else if (type == "density") {
        densityplot_match(x = x, t = t, w = w, sw = sw, ...)
      }

      devAskNewPage(ask = interactive)
    }
  }
  devAskNewPage(ask = FALSE)

 invisible(NULL)
}

qqplot_match <- function(x, t, w, sw, discrete.cutoff = 5, ...) {

  ord <- order(x)
  x_ord <- x[ord]
  t_ord <- t[ord]

  u <- unique(x_ord)

  #Need to interpolate larger group to be same size as smaller group

  #Unmatched sample
  sw_ord <- sw[ord]
  sw1 <- sw_ord[t_ord == 1]
  sw0 <- sw_ord[t_ord != 1]

  x1 <- x_ord[t_ord == 1][sw1 > 0]
  x0 <- x_ord[t_ord != 1][sw0 > 0]

  sw1 <- sw1[sw1 > 0]
  sw0 <- sw0[sw0 > 0]

  swn1 <- length(sw1)
  swn0 <- length(sw0)

  if (swn1 < swn0) {
    if (length(u) <= discrete.cutoff) {
      x0probs <- vapply(u, function(u_) wm(x0 == u_, sw0), numeric(1L))
      x0cumprobs <- c(0, .cumsum_prob(x0probs))
      x0 <- u[findInterval(.cumsum_prob(sw1), x0cumprobs, rightmost.closed = TRUE)]
    }
    else {
      x0 <- approx(.cumsum_prob(sw0), y = x0, xout = .cumsum_prob(sw1), rule = 2,
                   method = "constant", ties = "ordered")$y
    }
  }
  else if (swn1 > swn0) {
    if (length(u) <= discrete.cutoff) {
      x1probs <- vapply(u, function(u_) wm(x1 == u_, sw1), numeric(1L))
      x1cumprobs <- c(0, .cumsum_prob(x1probs))
      x1 <- u[findInterval(.cumsum_prob(sw0), x1cumprobs, rightmost.closed = TRUE)]
    }
    else {
      x1 <- approx(.cumsum_prob(sw1), y = x1, xout = .cumsum_prob(sw0), rule = 2,
                   method = "constant", ties = "ordered")$y
    }
  }

  if (length(u) <= discrete.cutoff) {
    md <- min(diff(u))
    x0 <- jitter(x0, amount = .1 * md)
    x1 <- jitter(x1, amount = .1 * md)
  }

  rr <- range(c(x0, x1))
  plot(x0, x1, xlab = "", ylab = "", xlim = rr, ylim = rr, axes = FALSE, ...)
  abline(a = 0, b = 1)
  abline(a = (rr[2]-rr[1])*0.1, b = 1, lty = 2)
  abline(a = -(rr[2]-rr[1])*0.1, b = 1, lty = 2)
  axis(2)
  box()

  #Matched sample
  w_ord <- w[ord]
  w1 <- w_ord[t_ord == 1]
  w0 <- w_ord[t_ord != 1]

  x1 <- x_ord[t_ord == 1][w1 > 0]
  x0 <- x_ord[t_ord != 1][w0 > 0]

  w1 <- w1[w1 > 0]
  w0 <- w0[w0 > 0]

  wn1 <- length(w1)
  wn0 <- length(w0)

  if (wn1 < wn0) {
    if (length(u) <= discrete.cutoff) {
      x0probs <- vapply(u, function(u_) wm(x0 == u_, w0), numeric(1L))
      x0cumprobs <- c(0, .cumsum_prob(x0probs))
      x0 <- u[findInterval(.cumsum_prob(w1), x0cumprobs, rightmost.closed = TRUE)]
    }
    else {
      x0 <- approx(.cumsum_prob(w0), y = x0, xout = .cumsum_prob(w1), rule = 2,
                   method = "constant", ties = "ordered")$y
    }
  }
  else if (wn1 > wn0) {
    if (length(u) <= discrete.cutoff) {
      x1probs <- vapply(u, function(u_) wm(x1 == u_, w1), numeric(1L))
      x1cumprobs <- c(0, .cumsum_prob(x1probs))
      x1 <- u[findInterval(.cumsum_prob(w0), x1cumprobs, rightmost.closed = TRUE)]
    }
    else {
      x1 <- approx(.cumsum_prob(w1), y = x1, xout = .cumsum_prob(w0), rule = 2,
                   method = "constant", ties = "ordered")$y
    }
  }

  if (length(u) <= discrete.cutoff) {
    md <- min(diff(u))
    x0 <- jitter(x0, amount = .1 * md)
    x1 <- jitter(x1, amount = .1 * md)
  }

  plot(x0, x1, xlab = "", ylab = "", xlim = rr, ylim = rr, axes = FALSE, ...)
  abline(a = 0, b = 1)
  abline(a = (rr[2]-rr[1])*0.1, b = 1, lty = 2)
  abline(a = -(rr[2]-rr[1])*0.1, b = 1, lty = 2)
  box()
}

ecdfplot_match <- function(x, t, w, sw, ...) {
  ord <- order(x)
  x.min <- x[ord][1]
  x.max <- x[ord][length(x)]
  x.range <- x.max - x.min

  #Unmatched samples
  plot(x = x, y = w, type= "n" , xlim = c(x.min - .02 * x.range, x.max + .02 * x.range),
       ylim = c(0, 1), axes = TRUE, ...)

  for (tr in 0:1) {
    in.tr <- t[ord] == tr
    ordt <- ord[in.tr]
    cswt <- c(0, .cumsum_prob(sw[ordt]), 1)
    xt <- c(x.min - .02 * x.range, x[ordt], x.max + .02 * x.range)

    lines(x = xt, y = cswt, type = "s", col = if (tr == 0) "grey60" else "black")
  }

  abline(h = 0:1)
  box()

  #Matched sample
  plot(x = x, y = w, type= "n" , xlim = c(x.min - .02 * x.range, x.max + .02 * x.range),
       ylim = c(0, 1), axes = FALSE, ...)
  for (tr in 0:1) {
    in.tr <- t[ord] == tr
    ordt <- ord[in.tr]
    cwt <- c(0, .cumsum_prob(w[ordt]), 1)
    xt <- c(x.min - .02 * x.range, x[ordt], x.max + .02 * x.range)

    lines(x = xt, y = cwt, type = "s", col = if (tr == 0) "grey60" else "black")

  }

  abline(h = 0:1)
  axis(1)
  box()
}

densityplot_match <- function(x, t, w, sw, ...) {

  if (length(unique(x)) == 2L) x <- factor(x)

  if (!is.factor(x)) {
    #Density plot for continuous variable
    small.tr <- (0:1)[which.min(c(sum(t==0), sum(t==1)))]
    x_small <- x[t==small.tr]

    x.min <- min(x)
    x.max <- max(x)
    #
    A <- list(...)

    bw <- A[["bw"]]
    if (is.null(bw)) A[["bw"]] <- bw.nrd0(x_small)
    else if (is.character(bw)) {
      bw <- tolower(bw)
      bw <- match_arg(bw, c("nrd0", "nrd", "ucv", "bcv", "sj", "sj-ste", "sj-dpi"))
      A[["bw"]] <- switch(bw, nrd0 = bw.nrd0(x_small), nrd = bw.nrd(x_small),
                          ucv = bw.ucv(x_small), bcv = bw.bcv(x_small), sj = ,
                          `sj-ste` = bw.SJ(x_small, method = "ste"),
                          `sj-dpi` = bw.SJ(x_small, method = "dpi"))
    }
    if (is.null(A[["cut"]])) A[["cut"]] <- 3

    d_unmatched <- do.call("rbind", lapply(0:1, function(tr) {
      cbind(as.data.frame(do.call("density", c(list(x[t==tr], weights = sw[t==tr],
                                                    from = x.min - A[["cut"]]*A[["bw"]],
                                                    to = x.max + A[["cut"]]*A[["bw"]]), A))[1:2]),
            t = tr)
    }))

    d_matched <- do.call("rbind", lapply(0:1, function(tr) {
      cbind(as.data.frame(do.call("density", c(list(x[t==tr], weights = w[t==tr],
                                                    from = x.min - A[["cut"]]*A[["bw"]],
                                                    to = x.max + A[["cut"]]*A[["bw"]]), A))[1:2]),
            t = tr)
    }))

    y.max <- max(d_unmatched$y, d_matched$y)

    #Unmatched samples
    plot(x = d_unmatched$x, y = d_unmatched$y, type = "n",
         xlim = c(x.min - A[["cut"]]*A[["bw"]], x.max + A[["cut"]]*A[["bw"]]),
         ylim = c(0, 1.1*y.max), axes = TRUE, ...)

    for (tr in 0:1) {
      in.tr <- d_unmatched$t == tr
      xt <- d_unmatched$x[in.tr]
      yt <- d_unmatched$y[in.tr]

      lines(x = xt, y = yt, type = "l", col = if (tr == 0) "grey60" else "black")
    }

    abline(h = 0)
    box()

    #Matched sample
    plot(x = d_matched$x, y = d_matched$y, type = "n",
         xlim = c(x.min - A[["cut"]]*A[["bw"]], x.max + A[["cut"]]*A[["bw"]]),
         ylim = c(0, 1.1*y.max), axes = FALSE, ...)

    for (tr in 0:1) {
      in.tr <- d_matched$t == tr
      xt <- d_matched$x[in.tr]
      yt <- d_matched$y[in.tr]

      lines(x = xt, y = yt, type = "l", col = if (tr == 0) "grey60" else "black")
    }

    abline(h = 0)
    axis(1)
    box()
  }
  else {
    #Bar plot for binary variable
    x_t_un <- lapply(sort(unique(t)), function(t_) {
      vapply(levels(x), function(i) {
      wm(x[t==t_] == i, sw[t==t_])
    }, numeric(1L))})

    x_t_m <- lapply(sort(unique(t)), function(t_) {
      vapply(levels(x), function(i) {
        wm(x[t==t_] == i, w[t==t_])
      }, numeric(1L))})

    ylim <- c(0, 1.1*max(unlist(x_t_un), unlist(x_t_m)))

    borders <- c("grey60", "black")
    for (i in seq_along(x_t_un)) {
      barplot(x_t_un[[i]], border = borders[i],
              col = if (i == 1) "white" else NA_character_,
              ylim = ylim, add = i != 1L)
    }

    abline(h = 0:1)
    box()

    for (i in seq_along(x_t_m)) {
      barplot(x_t_m[[i]], border = borders[i],
              col = if (i == 1) "white" else NA_character_,
              ylim = ylim, add = i != 1L, axes = FALSE)
    }

    abline(h = 0:1)
    box()
  }
}

hist.pscore <- function(x, xlab = "Propensity Score", freq = FALSE, ...){
  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  treat <- x$treat
  pscore <- x$distance[!is.na(x$distance)]
  s.weights <- if (is.null(x$s.weights)) rep(1, length(treat)) else x$s.weights
  weights <- x$weights * s.weights
  matched <- weights != 0
  q.cut <- x$q.cut

  minp <- min(pscore)
  maxp <- max(pscore)
  ratio <- x$call$ratio

  if (is.null(ratio)) ratio <- 1

  if (freq) {
    weights <- .make_sum_to_n(weights, by = treat)
    s.weights <- .make_sum_to_n(s.weights, by = treat)
  }
  else {
    weights <- .make_sum_to_1(weights, by = treat)
    s.weights <- .make_sum_to_1(s.weights, by = treat)
  }

  ylab <- if (freq) "Count" else "Proportion"

  par(mfrow = c(2,2))
  # breaks <- pretty(na.omit(pscore), 10)
  breaks <- seq(minp, maxp, length = 11)
  xlim <- range(breaks)

  for (n in c("Raw Treated", "Matched Treated", "Raw Control", "Matched Control")) {
    if (startsWith(n, "Raw")) w <- s.weights
    else w <- weights

    if (endsWith(n, "Treated")) t <- 1
    else t <- 0

    #Create histogram using weights
    #Manually assign density, which is used as height of the bars. The scaling
    #of the weights above determine whether they are "counts" or "proportions".
    #Regardless, set freq = FALSE in plot() to ensure density is used for bar
    #height rather than count.
    pm <- hist(pscore[treat==t], plot = FALSE, breaks = breaks)
    pm[["density"]] <- vapply(seq_len(length(pm$breaks) - 1), function(i) {
      sum(w[treat == t & pscore >= pm$breaks[i] & pscore < pm$breaks[i+1]])
    }, numeric(1L))
    plot(pm, xlim = xlim, xlab = xlab, main = n, ylab = ylab,
         freq = FALSE, col = "lightgray", ...)

    if (!startsWith(n, "Raw") && !is.null(q.cut)) abline(v = q.cut, lty=2)
  }

}

jitter.pscore <- function(x, interactive, pch = 1, ...){

  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  treat <- x$treat
  pscore <- x$distance
  s.weights <- if (is.null(x$s.weights)) rep(1, length(treat)) else x$s.weights
  weights <- x$weights * s.weights
  matched <- weights > 0
  q.cut <- x$q.cut
  jitp <- jitter(rep(1,length(treat)), factor=6)+(treat==1)*(weights==0)-(treat==0) - (weights==0)*(treat==0)
  cswt <- sqrt(s.weights)
  cwt <- sqrt(weights)
  minp <- min(pscore, na.rm = TRUE)
  maxp <- max(pscore, na.rm = TRUE)

  plot(pscore, xlim = c(minp - 0.05*(maxp-minp), maxp + 0.05*(maxp-minp)), ylim = c(-1.5,2.5),
       type = "n", ylab = "", xlab = "Propensity Score",
       axes = FALSE, main = "Distribution of Propensity Scores", ...)
  if (!is.null(q.cut)) abline(v = q.cut, col = "grey", lty = 1)

  #Matched treated
  points(pscore[treat==1 & matched], jitp[treat==1 & matched],
         pch = pch, cex = cwt[treat==1 & matched], ...)
  #Matched control
  points(pscore[treat==0 & matched], jitp[treat==0 & matched],
         pch = pch, cex = cwt[treat==0 & matched], ...)
  #Unmatched treated
  points(pscore[treat==1 & !matched], jitp[treat==1 & !matched],
         pch = pch, cex = cswt[treat==1 & matched],...)
  #Unmatched control
  points(pscore[treat==0 & !matched], jitp[treat==0 & !matched],
         pch = pch, cex = cswt[treat==0 & matched], ...)

  axis(1)

  center <- mean(par("usr")[1:2])
  text(center, 2.5, "Unmatched Treated Units", adj = .5)
  text(center, 1.5, "Matched Treated Units", adj = .5)
  text(center, 0.5, "Matched Control Units", adj = .5)
  text(center, -0.5, "Unmatched Control Units", adj = .5)
  box()

  if (interactive) {
    cat("To identify the units, use first mouse button; to stop, use second.\n")
    identify(pscore, jitp, names(treat), atpen = TRUE)
  }
}
