distance2rpart <- function(formula, data, ...) {
    if (requireNamespace("rpart", quietly = TRUE)) {
        res <- rpart::rpart(formula, data, ...)
        return(list(model = res, distance = predict(res)))
    } else {
        stop("rpart package is needed. Please install it.",
             call. = FALSE)
    }
}
       

