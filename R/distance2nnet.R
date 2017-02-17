distance2nnet <- function(formula, data, ...) {
    if (requireNamespace("nnet", quietly = TRUE)) {
        res <- nnet::nnet(formula, data, ...)
        return(list(model = res, distance = fitted(res)))
    } else {
        stop("nnet package is needed. Please install it.",
             call. = FALSE)
    }
}

