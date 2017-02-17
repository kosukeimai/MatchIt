distance2GAMlogit <- function(formula, data, ...) {
    if (requireNamespace("mgcv", quietly = TRUE)) {
        res <- mgcv::gam(formula, data, family=binomial(logit), ...)
        return(list(model = res, distance = fitted(res)))
    } else {
        stop("mgcv package is needed. Please install it.",
             call. = FALSE)
    }
}

distance2GAMprobit <- function(formula, data, ...) {
    if (requireNamespace("mgcv", quietly = TRUE)) {
        res <- mgcv::gam(formula, data, family=binomial(probit), ...)
        return(list(model = res, distance = fitted(res)))
    } else {
        stop("mgcv package is needed. Please install it.",
             call. = FALSE)
    }      
}

distance2GAMcloglog <- function(formula, data, ...) {
    if (requireNamespace("mgcv", quietly = TRUE)) {
        res <- mgcv::gam(formula, data, family=binomial(cloglog), ...)
        return(list(model = res, distance = fitted(res)))
    } else {
        stop("mgcv package is needed. Please install it.",
             call. = FALSE)
    }
}


distance2GAMlog <- function(formula, data, ...) {
    if (requireNamespace("mgcv", quietly = TRUE)) {
        res <- mgcv::gam(formula, data, family=binomial(log), ...)
        return(list(model = res, distance = fitted(res)))
    } else {
        stop("mgcv package is needed. Please install it.",
             call. = FALSE)
    }
}
        

distance2GAMcauchit <- function(formula, data, ...) {
    if (requireNamespace("mgcv", quietly = TRUE)) {
        res <- mgcv::gam(formula, data, family=binomial(cauchit), ...)
        return(list(model = res, distance = fitted(res)))
    } else {
        stop("mgcv package is needed. Please install it.",
             call. = FALSE)
    }
}
        

