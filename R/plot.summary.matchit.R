#' Generate a Love Plot of Standardized Mean Differences
#'
#' Generates a Love plot, which is a dot plot with variable names on the y-axis
#' and standardized mean differences on the x-axis. Each point represents the
#' standardized mean difference of the corresponding covariate in the matched
#' or unmatched sample. Love plots are a simple way to display covariate
#' balance before and after matching. The plots are generated using
#' [dotchart()] and [points()].
#'
#' @param x a `summary.matchit` object; the output of a call to
#' [summary.matchit()]. The `standardize` argument must be set to
#' `TRUE` (which is the default) in the call to `summary`.
#' @param abs `logical`; whether the standardized mean differences should
#' be displayed in absolute value (`TRUE`, default) or not `FALSE`.
#' @param var.order how the variables should be ordered. Allowable options
#' include `"data"`, ordering the variables as they appear in the
#' `summary` output; `"unmatched"`, ordered the variables based on
#' their standardized mean differences before matching; `"matched"`,
#' ordered the variables based on their standardized mean differences after
#' matching; and `"alphabetical"`, ordering the variables alphabetically.
#' Default is `"data"`. Abbreviations allowed.
#' @param threshold numeric values at which to place vertical lines indicating
#' a balance threshold. These can make it easier to see for which variables
#' balance has been achieved given a threshold. Multiple values can be supplied
#' to add multiple lines. When `abs = FALSE`, the lines will be displayed
#' on both sides of zero. The lines are drawn with `abline` with the
#' linetype (`lty`) argument corresponding to the order of the entered
#' variables (see options at [par()]). The default is `c(.1, .05)` for a
#' solid line (`lty = 1`) at .1 and a dashed line (`lty = 2`) at .05,
#' indicating acceptable and good balance, respectively. Enter a value as
#' `NA` to skip that value of `lty` (e.g., `c(NA, .05)` to have
#' only a dashed vertical line at .05).
#' @param position the position of the legend. Should be one of the allowed
#' keyword options supplied to `x` in [legend()] (e.g., `"right"`,
#' `"bottomright"`, etc.). Default is `"bottomright"`. Set to
#' `NULL` for no legend to be included. Note that the legend will cover up
#' points if you are not careful; setting `var.order` appropriately can
#' help in avoiding this.
#' @param \dots ignored.
#'
#' @return A plot is displayed, and `x` is invisibly returned.
#'
#' @details
#' For matching methods other than subclassification,
#' `plot.summary.matchit` uses `x$sum.all[,"Std. Mean Diff."]` and
#' `x$sum.matched[,"Std. Mean Diff."]` as the x-axis values. For
#' subclassification, in addition to points for the unadjusted and aggregate
#' subclass balance, numerals representing balance in individual subclasses are
#' plotted if `subclass = TRUE` in the call to `summary`. Aggregate
#' subclass standardized mean differences are taken from
#' `x$sum.across[,"Std. Mean Diff."]` and the subclass-specific mean
#' differences are taken from `x$sum.subclass`.
#'
#' @author Noah Greifer
#'
#' @seealso [summary.matchit()], [dotchart()]
#'
#' \pkgfun{cobalt}{love.plot} is a more flexible and sophisticated function to make
#' Love plots and is also natively compatible with `matchit` objects.
#'
#' @examples
#'
#' data("lalonde")
#' m.out <- matchit(treat ~ age + educ + married +
#'                    race + re74,
#'                  data = lalonde,
#'                  method = "nearest")
#' plot(summary(m.out, interactions = TRUE),
#'      var.order = "unmatched")
#'
#' s.out <- matchit(treat ~ age + educ + married +
#'                    race + nodegree + re74 + re75,
#'                  data = lalonde,
#'                  method = "subclass")
#' plot(summary(s.out, subclass = TRUE),
#'      var.order = "unmatched",
#'      abs = FALSE)
#'
#' @exportS3Method plot summary.matchit
plot.summary.matchit <- function(x,
                                 abs = TRUE,
                                 var.order = "data",
                                 threshold = c(.1, .05),
                                 position = "bottomright",
                                 ...) {

  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  sub <- inherits(x, "summary.matchit.subclass")
  matched <- sub || is_not_null(x[["sum.matched"]])
  un <- is_not_null(x[["sum.all"]])

  standard.sum <- {
    if (un) x[["sum.all"]]
    else if (sub) x[["sum.across"]]
    else x[["sum.matched"]]
  }

  if (!"Std. Mean Diff." %in% colnames(standard.sum)) {
    .err("not appropriate for unstandardized summary. Run `summary()` with the `standardize = TRUE` option, and then plot")
  }

  if (un) {
    sd.all <- x[["sum.all"]][, "Std. Mean Diff."]
  }

  if (matched) {
    sd.matched <- x[[if (sub) "sum.across" else "sum.matched"]][, "Std. Mean Diff."]
  }

  chk::chk_flag(abs)

  var.names <- rownames(standard.sum)

  chk::chk_string(var.order)
  var.order <- match_arg(var.order, c("data", "matched", "unmatched", "alphabetical"))

  if (!un && var.order == "unmatched") {
    .err('`var.order` cannot be "unmatched" if `un = TRUE` in the call to `summary()`')
  }

  if (!matched && var.order == "matched") {
    .err('`var.order` cannot be "matched" if `method = NULL` in the original call to `matchit()`')
  }

  if (abs) {
    if (un) sd.all <- abs(sd.all)
    if (matched) sd.matched <- abs(sd.matched)
    xlab <- "Absolute Standardized\nMean Difference"
  }
  else {
    xlab <- "Standardized Mean Difference"
  }

  ord <- switch(var.order,
                "data" = rev(seq_along(var.names)),
                "matched" = order(sd.matched),
                "unmatched" = order(sd.all),
                "alphabetical" = order(var.names, decreasing = TRUE))

  dotchart(if (un) sd.all[ord] else sd.matched[ord],
           labels = var.names[ord], xlab = xlab,
           bg = NA, color = NA, ...)
  abline(v = 0)

  if (sub && is_not_null(x$sum.subclass)) {
    for (i in seq_along(x$sum.subclass)) {
      sd.sub <- x$sum.subclass[[i]][, "Std. Mean Diff."]
      if (abs) sd.sub <- abs(sd.sub)
      points(x = sd.sub[ord], y = seq_along(sd.sub),
             pch = as.character(i), col = "gray60", cex = .6)
    }
  }

  if (un) {
    points(x = sd.all[ord], y = seq_along(sd.all),
           pch = 21, bg = "white", col = "black")
  }

  if (matched) {
    points(x = sd.matched[ord], y = seq_along(sd.matched),
           pch = 21, bg = "black", col = "black")
  }

  if (is_not_null(threshold)) {
    if (abs) {
      abline(v = threshold, lty = seq_along(threshold))
    }
    else {
      abline(v = threshold, lty = seq_along(threshold))
      abline(v = -threshold, lty = seq_along(threshold))
    }
  }

  if (sum(matched, un) > 1 && is_not_null(position)) {
    position <- match_arg(position, c("bottomright", "bottom", "bottomleft", "left",
                                      "topleft", "top", "topright", "right", "center"))
    legend(position, legend = c("All", "Matched"),
           pt.bg = c("white", "black"), pch = 21,
           inset = .015, xpd = TRUE)
  }

  invisible(x)
}
