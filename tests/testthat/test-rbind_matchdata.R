# Tests for rbind.matchdata() and rbind.getmatches().
#
# The job of these methods is to stack matched datasets that came from separate calls
# to matchit() without letting their subclasses collide or their differently named
# output columns drift apart, so that is what is checked here: the row and column
# arithmetic, the renaming, the subclass relabeling, and the errors that keep
# incompatible inputs from being silently combined.

data("lalonde", package = "MatchIt")

f_rb <- treat ~ age + educ + re74

#Two matched datasets from disjoint subsets of the same data, the use case the methods
#are written for. The black subset has more treated than control units, which warns.
expect_wrn(m_black <- matchit(f_rb, data = subset(lalonde, race == "black")),
           "fewer control units than treated units")

m_hispan <- matchit(f_rb, data = subset(lalonde, race == "hispan"))

md_black <- match_data(m_black)
md_hispan <- match_data(m_hispan)

test_that("rbind.matchdata: rows and columns add up", {
  all <- rbind(md_black, md_hispan)

  expect_s3_class(all, "matchdata")
  expect_s3_class(all, "data.frame")
  expect_identical(nrow(all), nrow(md_black) + nrow(md_hispan))
  expect_named(all, names(md_black))

  #The output columns carry the same values they had before stacking
  expect_equal(all$weights, c(md_black$weights, md_hispan$weights),
               ignore_attr = TRUE)
  expect_equal(all$distance, c(md_black$distance, md_hispan$distance),
               ignore_attr = TRUE)

  #And the attributes naming those columns survive
  expect_identical(attr(all, "distance"), "distance")
  expect_identical(attr(all, "weights"), "weights")
  expect_identical(attr(all, "subclass"), "subclass")
})

test_that("rbind.matchdata: subclasses are made unique across datasets", {
  all <- rbind(md_black, md_hispan)

  #Each input's levels are prefixed with its position, so no two inputs share a level
  expect_identical(levels(all$subclass),
                   c(paste0("1_", levels(md_black$subclass)),
                     paste0("2_", levels(md_hispan$subclass))))

  #Every unit keeps its own pair: the subclass sizes are the input sizes, relabeled
  expect_identical(unname(table(all$subclass)[paste0("1_", levels(md_black$subclass))]),
                   unname(table(md_black$subclass)))

  #No unit is left without a subclass
  expect_false(anyNA(all$subclass))
})

test_that("rbind.matchdata: differently named output columns are reconciled", {
  #The same information under a different name in one input
  md_hispan2 <- match_data(m_hispan, distance = "prop.score", weights = "w")

  all <- rbind(md_black, md_hispan2)

  #The first non-missing name wins, and both inputs' values end up under it
  expect_identical(attr(all, "distance"), "distance")
  expect_identical(attr(all, "weights"), "weights")
  expect_false(any(c("prop.score", "w") %in% names(all)))
  expect_false(anyNA(all$distance))
  expect_false(anyNA(all$weights))

  #Naming is decided by position, so the other order takes the other name
  all2 <- rbind(md_hispan2, md_black)
  expect_identical(attr(all2, "distance"), "prop.score")
  expect_identical(attr(all2, "weights"), "w")
  expect_true(all(c("prop.score", "w") %in% names(all2)))
})

test_that("rbind.matchdata: a column one input lacks is filled with NA", {
  #Coarsened exact matching produces no distance measure
  md_cem <- match_data(matchit(f_rb, data = subset(lalonde, race == "hispan"),
                               method = "cem"))
  expect_false(hasName(md_cem, "distance"))

  all <- rbind(md_black, md_cem)

  expect_true(hasName(all, "distance"))
  expect_identical(sum(is.na(all$distance)), nrow(md_cem))
  expect_false(anyNA(all$distance[seq_len(nrow(md_black))]))
})

test_that("rbind.matchdata: column order does not have to agree", {
  #The columns of the second input are reordered to match the first
  reordered <- md_hispan[rev(names(md_hispan))]
  class(reordered) <- class(md_hispan)
  attributes(reordered)[c("distance", "weights", "subclass")] <-
    attributes(md_hispan)[c("distance", "weights", "subclass")]

  all <- rbind(md_black, reordered)

  expect_named(all, names(md_black))
  expect_identical(nrow(all), nrow(md_black) + nrow(md_hispan))
})

test_that("rbind.matchdata: more than two datasets can be stacked", {
  m_white <- matchit(f_rb, data = subset(lalonde, race == "white"))

  md_white <- match_data(m_white)

  all <- rbind(md_black, md_hispan, md_white)

  expect_identical(nrow(all), nrow(md_black) + nrow(md_hispan) + nrow(md_white))
  expect_identical(nlevels(all$subclass),
                   nlevels(md_black$subclass) + nlevels(md_hispan$subclass) +
                     nlevels(md_white$subclass))

  #Levels are prefixed by position, so the third input's are prefixed with 3
  expect_true(all(startsWith(levels(all$subclass)[nlevels(all$subclass)], "3_")))
})

test_that("rbind.matchdata: rejects inputs that should not be combined", {
  g_black <- get_matches(m_black)

  #matchdata and getmatches are not interchangeable
  expect_err(rbind(md_black, g_black),
             "supplied objects must be all <matchdata> objects or all <getmatches> objects")

  #Different columns mean different datasets
  expect_err(rbind(md_black, md_black[-2L]),
             "the `match_data()` inputs must come from the same dataset")

  #Nothing to dispatch on
  expect_err(rbind.matchdata(lalonde),
             "a <matchdata> or <getmatches> object must be supplied")
})

test_that("rbind.getmatches: works the same way and keeps the id column", {
  g_black <- get_matches(m_black)
  g_hispan <- get_matches(m_hispan)

  all <- rbind(g_black, g_hispan)

  expect_s3_class(all, "getmatches")
  expect_identical(nrow(all), nrow(g_black) + nrow(g_hispan))
  expect_identical(attr(all, "id"), "id")
  expect_identical(attr(all, "subclass"), "subclass")

  #Subclasses are relabeled here too, so pairs from different fits stay distinct
  expect_identical(nlevels(all$subclass),
                   nlevels(g_black$subclass) + nlevels(g_hispan$subclass))

  #A renamed id is reconciled like the other output columns
  g_hispan2 <- get_matches(m_hispan, id = "unit")
  all2 <- rbind(g_black, g_hispan2)
  expect_identical(attr(all2, "id"), "id")
  expect_false(hasName(all2, "unit"))
  expect_false(anyNA(all2$id))

  #And the error for mixing types points the other way
  expect_err(rbind(g_black, md_black),
             "supplied objects must be all <matchdata> objects or all <getmatches> objects")
})

test_that("rbind.matchdata: the output can be stacked again", {
  #`rbind()`'s result is itself a matchdata object, so it dispatches to the same method
  m_white <- matchit(f_rb, data = subset(lalonde, race == "white"))

  once <- rbind(md_black, md_hispan)
  twice <- rbind(once, match_data(m_white))

  expect_s3_class(twice, "matchdata")
  expect_identical(nrow(twice), nrow(once) + sum(m_white$weights > 0))

  #Subclasses from the already-stacked input keep their old prefix under a new one
  expect_true(any(startsWith(levels(twice$subclass), "1_1_")))
  expect_false(anyNA(twice$subclass))
})

test_that("rbind.matchdata: a single input is returned essentially unchanged", {
  one <- rbind(md_black)

  expect_s3_class(one, "matchdata")
  expect_identical(nrow(one), nrow(md_black))
  expect_named(one, names(md_black))

  #Only the subclass labels change, gaining the position prefix
  expect_identical(levels(one$subclass), paste0("1_", levels(md_black$subclass)))
  expect_identical(as.data.frame(one)[setdiff(names(one), "subclass")],
                   as.data.frame(md_black)[setdiff(names(md_black), "subclass")])
})

test_that("rbind.matchdata: emits no conditions", {
  expect_silent(rbind(md_black, md_hispan))
  expect_silent(rbind(get_matches(m_black), get_matches(m_hispan)))
})
