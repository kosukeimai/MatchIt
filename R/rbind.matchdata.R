#' Append matched datasets together
#'
#' These functions are [rbind()] methods for objects resulting from calls to
#' [match.data()] and [get_matches()]. They function nearly identically to
#' `rbind.data.frame()`; see Details for how they differ.
#'
#' @aliases rbind.matchdata rbind.getmatches
#'
#' @param \dots Two or more `matchdata` or `getmatches` objects the
#' output of calls to [match.data()] and [get_matches()], respectively.
#' Supplied objects must either be all `matchdata` objects or all
#' `getmatches` objects.
#' @param deparse.level Passed to [rbind()].
#'
#' @return An object of the same class as those supplied to it (i.e., a
#' `matchdata` object if `matchdata` objects are supplied and a
#' `getmatches` object if `getmatches` objects are supplied).
#' [rbind()] is called on the objects after adjusting the variables so that the
#' appropriate method will be dispatched corresponding to the class of the
#' original data object.
#'
#' @details
#' `rbind()` appends two or more datasets row-wise. This can be useful
#' when matching was performed separately on subsets of the original data and
#' they are to be combined into a single dataset for effect estimation. Using
#' the regular `data.frame` method for `rbind()` would pose a
#' problem, however; the `subclass` variable would have repeated names
#' across different datasets, even though units only belong to the subclasses
#' in their respective datasets. `rbind.matchdata()` renames the
#' subclasses so that the correct subclass membership is maintained.
#'
#' The supplied matched datasets must be generated from the same original
#' dataset, that is, having the same variables in it. The added components
#' (e.g., weights, subclass) can be named differently in different datasets but
#' will be changed to have the same name in the output.
#'
#' `rbind.getmatches()` and `rbind.matchdata()` are identical.
#'
#' @author Noah Greifer
#' @seealso [match.data()], [rbind()]
#'
#' See `vignettes("estimating-effects")` for details on using
#' `rbind()` for effect estimation after subsetting the data.
#'
#' @examples
#'
#' data("lalonde")
#'
#' # Matching based on race subsets
#' m.out_b <- matchit(treat ~ age + educ + married +
#'                     nodegree + re74 + re75,
#'                   data = subset(lalonde, race == "black"))
#' md_b <- match.data(m.out_b)
#'
#' m.out_h <- matchit(treat ~ age + educ + married +
#'                     nodegree + re74 + re75,
#'                   data = subset(lalonde, race == "hispan"))
#' md_h <- match.data(m.out_h)
#'
#' m.out_w <- matchit(treat ~ age + educ + married +
#'                     nodegree + re74 + re75,
#'                   data = subset(lalonde, race == "white"))
#' md_w <- match.data(m.out_w)
#'
#' #Bind the datasets together
#' md_all <- rbind(md_b, md_h, md_w)
#'
#' #Subclass conflicts are avoided
#' levels(md_all$subclass)
#'
#' @exportS3Method rbind matchdata
rbind.matchdata <- function(..., deparse.level = 1) {

  allargs <- list(...)
  allargs <- allargs[lengths(allargs) > 0L]
  if (is.null(names(allargs))) {
    md_list <- allargs
    allargs <- list()
  }
  else {
    md_list <- allargs[names(allargs) == ""]
    allargs[names(allargs) == ""] <- NULL
  }
  allargs$deparse.level <- deparse.level

  type <- intersect(c("matchdata", "getmatches"), unlist(lapply(md_list, class)))
  if (length(type) == 0) stop("A matchdata or getmatches object must be supplied.", call. = FALSE)
  else if (length(type) == 2) stop("Supplied objects must be all matchdata objects or all getmatches objects.", call. = FALSE)

  attrs <- c("distance", "weights", "subclass", "id")
  attr_list <- setNames(vector("list", length(attrs)), attrs)
  key_attrs <- setNames(rep(NA_character_, length(attrs)), attrs)

  for (i in attrs) {
    attr_list[[i]] <- unlist(lapply(md_list, function(m) {
      a <- attr(m, i)
      if (length(a) == 0) NA_character_ else a
    }))
    if (all(is.na(attr_list[[i]]))) attr_list[[i]] <- NULL
    else key_attrs[i] <- attr_list[[i]][which(!is.na(attr_list[[i]]))[1]]
  }
  attrs <- names(attr_list)
  key_attrs <- key_attrs[attrs]

  #Check if all non-attr columns are the same across datasets
  other_col_list <- lapply(seq_along(md_list), function(d) {
    setdiff(names(md_list[[d]]), unlist(lapply(attr_list, `[`, d)))
  })
  for (d in seq_along(md_list)[-1]) {
    if (length(other_col_list[[d]]) != length(other_col_list[[1]]) || !all(other_col_list[[d]] %in% other_col_list[[1]])) {
      stop(paste("The", switch(type, "matchdata" = "match.data()", "get_matches()"), "inputs must come from the same dataset."), call. = FALSE)
    }
  }

  for (d in seq_along(md_list)) {
    for (i in attrs) {
      #Rename columns of each attribute the same across datasets
      if (is.null(attr(md_list[[d]], i))) md_list[[d]] <- setNames(cbind(md_list[[d]], NA), c(names(md_list[[d]]), key_attrs[i]))
      else names(md_list[[d]])[names(md_list[[d]]) == attr_list[[i]][d]] <- key_attrs[i]

      #Give subclasses unique values across datasets
      if (i == "subclass") {
        if (all(is.na(md_list[[d]][[key_attrs[i]]]))) md_list[[d]][[key_attrs[i]]] <- factor(md_list[[d]][[key_attrs[i]]], levels = NA)
        else levels(md_list[[d]][[key_attrs[i]]]) <- paste(d, levels(md_list[[d]][[key_attrs[i]]]), sep = "_")
      }
    }

    #Put all columns in the same order
    if (d > 1) {
      md_list[[d]] <- md_list[[d]][names(md_list[[1]])]
    }
    class(md_list[[d]]) <- class(md_list[[d]])[class(md_list[[d]]) != type]
  }

  out <- do.call("rbind", c(md_list, allargs))

  for (i in attrs) {
    attr(out, i) <- unname(key_attrs[i])
  }

  class(out) <- c(type, class(out))

  out
}

#' @exportS3Method rbind getmatches
#' @rdname rbind.matchdata
rbind.getmatches <- rbind.matchdata
