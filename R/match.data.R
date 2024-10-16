#' Construct a matched dataset from a `matchit` object
#' @name match.data
#' @aliases match.data get_matches
#'
#' @description
#' `match.data()` and `get_matches()` create a data frame with
#' additional variables for the distance measure, matching weights, and
#' subclasses after matching. This dataset can be used to estimate treatment
#' effects after matching or subclassification. `get_matches()` is most
#' useful after matching with replacement; otherwise, `match.data()` is
#' more flexible. See Details below for the difference between them.
#'
#' @param object a `matchit` object; the output of a call to [matchit()].
#' @param group which group should comprise the matched dataset: `"all"`
#' for all units, `"treated"` for just treated units, or `"control"`
#' for just control units. Default is `"all"`.
#' @param distance a string containing the name that should be given to the
#' variable containing the distance measure in the data frame output. Default
#' is `"distance"`, but `"prop.score"` or similar might be a good
#' alternative if propensity scores were used in matching. Ignored if a
#' distance measure was not supplied or estimated in the call to
#' `matchit()`.
#' @param weights a string containing the name that should be given to the
#' variable containing the matching weights in the data frame output. Default
#' is `"weights"`.
#' @param subclass a string containing the name that should be given to the
#' variable containing the subclasses or matched pair membership in the data
#' frame output. Default is `"subclass"`.
#' @param id a string containing the name that should be given to the variable
#' containing the unit IDs in the data frame output. Default is `"id"`.
#' Only used with `get_matches()`; for `match.data()`, the units IDs
#' are stored in the row names of the returned data frame.
#' @param data a data frame containing the original dataset to which the
#' computed output variables (`distance`, `weights`, and/or
#' `subclass`) should be appended. If empty, `match.data()` and
#' `get_matches()` will attempt to find the dataset using the environment
#' of the `matchit` object, which can be unreliable; see Notes.
#' @param include.s.weights `logical`; whether to multiply the estimated
#' weights by the sampling weights supplied to `matchit()`, if any.
#' Default is `TRUE`. If `FALSE`, the weights in the
#' `match.data()` or `get_matches()` output should be multiplied by
#' the sampling weights before being supplied to the function estimating the
#' treatment effect in the matched data.
#' @param drop.unmatched `logical`; whether the returned data frame should
#' contain all units (`FALSE`) or only units that were matched (i.e., have
#' a matching weight greater than zero) (`TRUE`). Default is `TRUE`
#' to drop unmatched units.
#'
#' @details
#' `match.data()` creates a dataset with one row per unit. It will be
#' identical to the dataset supplied except that several new columns will be
#' added containing information related to the matching. When
#' `drop.unmatched = TRUE`, the default, units with weights of zero, which
#' are those units that were discarded by common support or the caliper or were
#' simply not matched, will be dropped from the dataset, leaving only the
#' subset of matched units. The idea is for the output of `match.data()`
#' to be used as the dataset input in calls to `glm()` or similar to
#' estimate treatment effects in the matched sample. It is important to include
#' the weights in the estimation of the effect and its standard error. The
#' subclass column, when created, contains pair or subclass membership and
#' should be used to estimate the effect and its standard error. Subclasses
#' will only be included if there is a `subclass` component in the
#' `matchit` object, which does not occur with matching with replacement,
#' in which case `get_matches()` should be used. See
#' `vignette("estimating-effects")` for information on how to use
#' `match.data()` output to estimate effects.
#'
#' `get_matches()` is similar to `match.data()`; the primary
#' difference occurs when matching is performed with replacement, i.e., when
#' units do not belong to a single matched pair. In this case, the output of
#' `get_matches()` will be a dataset that contains one row per unit for
#' each pair they are a part of. For example, if matching was performed with
#' replacement and a control unit was matched to two treated units, that
#' control unit will have two rows in the output dataset, one for each pair it
#' is a part of. Weights are computed for each row, and, for control units, are equal to the
#' inverse of the number of control units in each control unit's subclass; treated units get a weight of 1.
#' Unmatched units are dropped. An additional column with unit IDs will be
#' created (named using the `id` argument) to identify when the same unit
#' is present in multiple rows. This dataset structure allows for the inclusion
#' of both subclass membership and repeated use of units, unlike the output of
#' `match.data()`, which lacks subclass membership when matching is done
#' with replacement. A `match.matrix` component of the `matchit`
#' object must be present to use `get_matches()`; in some forms of
#' matching, it is absent, in which case `match.data()` should be used
#' instead. See `vignette("estimating-effects")` for information on how to
#' use `get_matches()` output to estimate effects after matching with
#' replacement.
#'
#' @return
#' A data frame containing the data supplied in the `data` argument or in the
#' original call to `matchit()` with the computed
#' output variables appended as additional columns, named according the
#' arguments above. For `match.data()`, the `group` and
#' `drop.unmatched` arguments control whether only subsets of the data are
#' returned. See Details above for how `match.data()` and
#' `get_matches()` differ. Note that `get_matches` sorts the data by
#' subclass and treatment status, unlike `match.data()`, which uses the
#' order of the data.
#'
#' The returned data frame will contain the variables in the original data set
#' or dataset supplied to `data` and the following columns:
#'
#' \item{distance}{The propensity score, if estimated or supplied to the
#' `distance` argument in `matchit()` as a vector.}
#' \item{weights}{The computed matching weights. These must be used in effect
#' estimation to correctly incorporate the matching.}
#' \item{subclass}{Matching
#' strata membership. Units with the same value are in the same stratum.}
#' \item{id}{The ID of each unit, corresponding to the row names in the
#' original data or dataset supplied to `data`. Only included in
#' `get_matches` output. This column can be used to identify which rows
#' belong to the same unit since the same unit may appear multiple times if
#' reused in matching with replacement.}
#'
#' These columns will take on the name supplied to the corresponding arguments
#' in the call to `match.data()` or `get_matches()`. See Examples for
#' an example of rename the `distance` column to `"prop.score"`.
#'
#' If `data` or the original dataset supplied to `matchit()` was a
#' `data.table` or `tbl`, the `match.data()` output will have
#' the same class, but the `get_matches()` output will always be a base R
#' `data.frame`.
#'
#' In addition to their base class (e.g., `data.frame` or `tbl`),
#' returned objects have the class `matchdata` or `getmatches`. This
#' class is important when using [`rbind()`][rbind.matchdata] to
#' append matched datasets.
#'
#' @note The most common way to use `match.data()` and
#' `get_matches()` is by supplying just the `matchit` object, e.g.,
#' as `match.data(m.out)`. A data set will first be searched in the
#' environment of the `matchit` formula, then in the calling environment
#' of `match.data()` or `get_matches()`, and finally in the
#' `model` component of the `matchit` object if a propensity score
#' was estimated.
#'
#' When called from an environment different from the one in which
#' `matchit()` was originally called and a propensity score was not
#' estimated (or was but with `discard` not `"none"` and
#' `reestimate = TRUE`), this syntax may not work because the original
#' dataset used to construct the matched dataset will not be found. This can
#' occur when `matchit()` was run within an [lapply()] or
#' `purrr::map()` call. The solution, which is recommended in all cases,
#' is simply to supply the original dataset to the `data` argument of
#' `match.data()`, e.g., as `match.data(m.out, data = original_data)`, as demonstrated in the Examples.
#'
#' @seealso
#'
#' [matchit()]; [rbind.matchdata()]
#'
#' `vignette("estimating-effects")` for uses of `match.data()` and
#' `get_matches()` in estimating treatment effects.
#'
#' @examples
#'
#' data("lalonde")
#'
#' # 4:1 matching w/replacement
#' m.out1 <- matchit(treat ~ age + educ + married +
#'                     race + nodegree + re74 + re75,
#'                   data = lalonde, replace = TRUE,
#'                   caliper = .05, ratio = 4)
#'
#' m.data1 <- match.data(m.out1, data = lalonde,
#'                       distance = "prop.score")
#' dim(m.data1) #one row per matched unit
#' head(m.data1, 10)
#'
#' g.matches1 <- get_matches(m.out1, data = lalonde,
#'                           distance = "prop.score")
#' dim(g.matches1) #multiple rows per matched unit
#' head(g.matches1, 10)
#'

#' @export
match.data <- function(object, group = "all", distance = "distance", weights = "weights", subclass = "subclass",
                       data = NULL, include.s.weights = TRUE, drop.unmatched = TRUE) {

  chk::chk_is(object, "matchit")

  if (is.null(data)) {
    env <- environment(object$formula)
    data <- try(eval(object$call$data, envir = env), silent = TRUE)
    if (length(data) == 0 || inherits(data, "try-error") ||
        length(dim(data)) != 2 || nrow(data) != length(object[["treat"]])) {
      env <- parent.frame()
      data <- try(eval(object$call$data, envir = env), silent = TRUE)
      if (length(data) == 0 || inherits(data, "try-error") ||
          length(dim(data)) != 2 || nrow(data) != length(object[["treat"]])) {
        data <- object[["model"]][["data"]]
        if (length(data) == 0 || nrow(data) != length(object[["treat"]])) {
          .err("a valid dataset could not be found. Please supply an argument to `data` containing the original dataset used in the matching")
        }
      }
    }
  }

  if (!is.data.frame(data)) {
    if (!is.matrix(data)) {
      .err("`data` must be a data frame")
    }
    data <- as.data.frame.matrix(data)
  }

  if (nrow(data) != length(object$treat)) {
    .err("`data` must have as many rows as there were units in the original call to `matchit()`")
  }

  if (!is.null(object$distance)) {
    chk::chk_not_null(distance)
    chk::chk_string(distance)
    if (distance %in% names(data)) {
      .err(sprintf("%s is already the name of a variable in the data. Please choose another name for distance using the `distance` argument",
                   add_quotes(distance)))
    }
    data[[distance]] <- object$distance
  }

  if (!is.null(object$weights)) {
    chk::chk_not_null(weights)
    chk::chk_string(weights)
    if (weights %in% names(data)) {
      .err(sprintf("%s is already the name of a variable in the data. Please choose another name for weights using the `weights` argument",
                   add_quotes(weights)))
    }
    data[[weights]] <- object$weights

    if (!is.null(object$s.weights) && include.s.weights) {
      data[[weights]] <- data[[weights]] * object$s.weights
    }
  }

  if (!is.null(object$subclass)) {
    chk::chk_not_null(subclass)
    chk::chk_string(subclass)
    if (subclass %in% names(data)) {
      .err(sprintf("%s is already the name of a variable in the data. Please choose another name for subclass using the `subclass` argument",
                   add_quotes(subclass)))
    }
    data[[subclass]] <- object$subclass
  }

  treat <- object$treat

  if (drop.unmatched && !is.null(object$weights)) {
    data <- data[object$weights > 0,,drop = FALSE]
    treat <- treat[object$weights > 0]
  }

  group <- match_arg(group, c("all", "treated", "control"))
  if (group == "treated") data <- data[treat == 1,,drop = FALSE]
  else if (group == "control") data <- data[treat == 0,,drop = FALSE]

  if (!is.null(object$distance)) attr(data, "distance") <- distance
  if (!is.null(object$weights)) attr(data, "weights") <- weights
  if (!is.null(object$subclass)) attr(data, "subclass") <- subclass

  class(data) <- c("matchdata", class(data))

  data
}

#' @export
#' @rdname match.data
get_matches <- function(object, distance = "distance", weights = "weights", subclass = "subclass",
                        id = "id", data = NULL, include.s.weights = TRUE) {

  chk::chk_is(object, "matchit")

  if (is.null(object$match.matrix)) {
    .err("a match.matrix component must be present in the matchit object, which does not occur with all types of matching. Use `match.data()` instead")
  }

  #Get initial data using match.data; note weights and subclass will be removed,
  #including them here just checks their names don't clash
  m.data <- match.data(object, group = "all", distance = distance,
                       weights = weights, subclass = subclass, data = data,
                       include.s.weights = FALSE, drop.unmatched = TRUE)

  chk::chk_not_null(id)
  chk::chk_string(id)

  if (id %in% names(m.data)) {
    .err(sprintf("%s is already the name of a variable in the data. Please choose another name for id using the `id` argument",
                 add_quotes(id)))
  }

  m.data[[id]] <- names(object$treat)[object$weights > 0]

  for (i in c(weights, subclass)) {
    if (i %in% names(m.data)) m.data[[i]] <- NULL
  }

  mm <- object$match.matrix
  mm <- mm[!is.na(mm[,1]),,drop = FALSE]
  tmm <- t(mm)

  num.matches <- rowSums(!is.na(mm))

  matched <- as.data.frame(matrix(NA_character_, nrow = nrow(mm) + sum(!is.na(mm)), ncol = 3))
  names(matched) <- c(id, subclass, weights)

  matched[[id]] <- c(as.vector(tmm[!is.na(tmm)]), rownames(mm))
  matched[[subclass]] <- c(as.vector(col(tmm)[!is.na(tmm)]), seq_len(nrow(mm)))
  matched[[weights]] <- c(1/num.matches[matched[[subclass]][seq_len(sum(!is.na(mm)))]], rep(1, nrow(mm)))

  if (!is.null(object$s.weights) && include.s.weights) {
    matched[[weights]] <- matched[[weights]] * object$s.weights[matched[[id]]]
  }

  out <- merge(matched, m.data, by = id, all.x = TRUE, sort = FALSE)

  out <- out[order(out[[subclass]], object$treat[out[[id]]], method = "radix", decreasing = c(FALSE, TRUE)),]
  rownames(out) <- NULL

  out[[subclass]] <- factor(out[[subclass]], labels = seq_len(nrow(mm)))

  if (!is.null(object$distance)) attr(out, "distance") <- distance
  attr(out, "weights") <- weights
  attr(out, "subclass") <- subclass
  attr(out, "id") <- id

  class(out) <- c("getmatches", class(out))

  out
}
