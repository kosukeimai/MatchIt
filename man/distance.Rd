\name{distance}
\alias{distance}

\title{Propensity scores and other distance measures}

\description{
Several matching methods require or can involve the distance between treated and control units. Options include the Mahalanobis distance, propensity score distance, or distance between user-supplied values. Propensity scores are also used for common support via the \code{discard} options and for defined calipers. This page documents the options that can be supplied to the \code{distance} argument to \code{\link{matchit}}.

There are two ways to specify the \code{distance} argument: 1) as the string \code{"mahalanobis"}, 2) as a string containing the name of a method for estimating propensity scores, or 3) as a vector of values whose pairwise differences define the distance between units.

When \code{distance} is specified as one of the allowed strings (described below) other than \code{"mahalanobis"}, a propensity score is estimated using the variables in \code{formula} and the method corresponding to the given argument. This propensity score can be used to compute the distance between units as the absolute difference between the propensity scores of pairs of units. In this respect, the propensity score is more like a "position" measure than a distance measure, since it is the pairwise difference that form the distance rather than the propensity scores themselves. Still, this naming convention is used to reflect their primary purpose without committing to the status of the estimated values as propensity scores, since transformations of the scores are allowed and user-supplied values that are not propensity scores can also be supplied (detailed below). Propensity scores can also be used to create calipers and common support restrictions, whether or not they are used in the actual distance measure used in the matching, if any.

In addition to the \code{distance} argument, two other arguments can be specified that relate to the estimation and manipulation of the propensity scores. The \code{link} argument allows for different links to be used in models that require them such as generalized linear models, for which the logit and probit links are allowed, among others. In addition to specifying the link, the \code{link} argument can be used to specify whether the propensity score or the linearized version of the propensity score should be used; by specifying \code{link = "linear.{link}"}, the linearized version will be used.

The \code{distance.options} argument can also be specified, which should be a list of values passed to the propensity score-estimating function, for example, to choose specific options or tuning parameters for the estimation method. If \code{formula}, \code{data}, or \code{verbose} are not supplied to \code{distance.options}, the corresponding arguments from \code{matchit} will be automatically supplied. See the Examples for demonstrations of the uses of \code{link} and \code{distance.options}. When \code{s.weights} is supplied in the call to \code{matchit}, it will automatically be passed to the propensity score-estimating function as the \code{weights} argument unless otherwise described below.
}

\section{Allowable options}{
Below are the allowed options for \code{distance}:

\describe{
  \item{\code{"glm"}}{
The propensity scores are estimated using a generalized linear model (e.g., logistic regression). The \code{formula} supplied to \code{matchit} is passed directly to \code{\link{glm}}, and \code{\link{predict.glm}} is used to compute the propensity scores. The \code{link} argument can be specified as a link function supplied to \code{\link{binomial}}, e.g., \code{"logit"}, which is the default. When \code{link} is prepended by \code{"linear."}, the linear predictor is used instead of the predicted probabilities. \code{distance = "glm"} with \code{link = "logit"} (logistic regression) is the default in \code{matchit}.
}
  \item{\code{"gam"}}{
The propensity scores are estimated using a generalized additive model. The \code{formula} supplied to \code{matchit} is passed directly to \pkgfun{mgcv}{gam} in the \pkg{mgcv} package, and \pkgfun{mgcv}{predict.gam} is used to compute the propensity scores. The \code{link} argument can be specified as a link function supplied to \code{\link{binomial}}, e.g., \code{"logit"}, which is the default. When \code{link} is prepended by \code{"linear."}, the linear predictor is used instead of the predicted probabilities. Note that unless the smoothing functions \pkgfun{mgcv}{s}, \pkgfun{mgcv}{te}, \pkgfun{mgcv}{ti}, or \pkgfun{mgcv}{t2} are used in \code{formula}, a generalized additive model is identical to a generalized linear model and will estimate the same propensity scores as \code{glm}. See the documentation for \pkgfun{mgcv}{gam}, \pkgfun{mgcv}{formula.gam}, and \\pkgfun{mgcv}{gam.models} for more information on how to specify these models. Also note that the formula returned in the \code{matchit} output object will be a simplified version of the supplied formula with smoothing terms removed (but all named variables present).
}
  \item{\code{"rpart"}}{
The propensity scores are estimated using a classification tree. The \code{formula} supplied to \code{matchit} is passed directly to \pkgfun{rpart}{rpart} in the \pkg{rpart} package, and \pkgfun{rpart}{predict.rpart} is used to compute the propensity scores. The \code{link} argument is ignored, and predicted probabilities are always returned as the distance measure.
}
  \item{\code{"randomforest"}}{
The propensity scores are estimated using a random forest. The \code{formula} supplied to \code{matchit} is passed directly to \pkgfun{randomForest}{randomForest} in the \pkg{randomForest} package, and \pkgfun{randomForest}{predict.randomForest} is used to compute the propensity scores. The \code{link} argument is ignored, and predicted probabilities are always returned as the distance measure. When \code{s.weights} is supplied to \code{matchit}, it will not be passed to \code{randomForest} because \code{randomForest} does not accept weights.
}
  \item{\code{"nnet"}}{
The propensity scores are estimated using a single-hidden-layer neural network. The \code{formula} supplied to \code{matchit} is passed directly to \pkgfun{nnet}{nnet} in the \pkg{nnet} package, and \code{\link{fitted}} is used to compute the propensity scores. The \code{link} argument is ignored, and predicted probabilities are always returned as the distance measure. An argument to \code{size} must be supplied to \code{distance.options} when using \code{method = "nnet"}.
}
  \item{\code{"cbps"}}{
The propensity scores are estimated using the covariate balancing propensity score (CBPS) algorithm, which is a form of logistic regression where balance constraints are incorporated to a generalized method of moments estimation of of the model coefficients. The \code{formula} supplied to \code{matchit} is passed directly to \pkgfun{CBPS}{CBPS} in the \pkg{CBPS} package, and \code{\link{fitted}} is used to compute the propensity scores. The \code{link} argument can be specified as \code{"linear"} to use the linear predictor instead of the predicted probabilities. No other links are allowed. The \code{estimand} argument supplied to \code{matchit} will be used to select the appropriate estimand for use in defining the balance constraints, so no argument needs to be supplied to \code{ATT} in \code{CBPS}.
}
  \item{\code{"bart"}}{
The propensity scores are estimated using Bayesian additive regression trees (BART). The \code{formula} supplied to \code{matchit} is passed directly to \pkgfun2{dbarts}{bart}{bart2} in the \pkg{dbarts} package, and \pkgfun2{dbarts}{bart}{fitted} is used to compute the propensity scores. The \code{link} argument can be specified as \code{"linear"} to use the linear predictor instead of the predicted probabilities. When \code{s.weights} is supplied to \code{matchit}, it will not be passed to \code{bart2} because the \code{weights} argument in \code{bart2} does not correspond to sampling weights.
}
  \item{\code{"mahalanobis"}}{
No propensity scores are estimated. Rather than using the propensity score difference as the distance between units, the Mahalanobis distance is used instead. See \code{\link{mahalanobis}} for details on how it is computed. The Mahalanobis distance is always computed using all the variables in \code{formula}. With this specification, calipers and common support restrictions cannot be used and the \code{distance} component of the output object will be empty because no propensity scores are estimated. The \code{link} and \code{distance.options} arguments are ignored. See individual methods pages for whether the Mahalanobis distance is allowed and how it is used. Sometimes this setting is just a placeholder to indicate that no propensity score is to be estimated (e.g., with \code{method = "genetic"}). To perform Mahalanobis distance matching \emph{and} estimate propensity scores to be used for a purpose other than matching, the \code{mahvars} argument should be used along with a different specification to \code{distance}. See the individual matching method pages for details on how to use \code{mahvars}.
}
}

\code{distance} can also be supplied as a numeric vector whose values will be taken to function like propensity scores; their pairwise difference will define the distance between units. This might be beneficial to supply propensity scores computed outside \code{matchit} or to resupply \code{matchit} with propensity scores estimated before without having to recompute them. When \code{distance} is a supplied as a numeric vector, \code{link} and \code{distance.options} are ignored.
}
\section{Outputs}{
When specifying an argument to \code{distance} that estimates a propensity score, the output of the function called to estimate the propensity score (e.g., the \code{glm} object when \code{distance = "glm"}) will be included in the \code{matchit} output object in the \code{model} component. When \code{distance} is anything other than \code{"mahalanobis"}, the estimated or supplied distance measures will be included in the \code{matchit} output object in the \code{distance} component.
}
\note{
In versions of \pkg{MatchIt} prior to 4.0.0, \code{distance} was specified in a slightly different way. When specifying arguments using the old syntax, they will automatically be converted to the corresponding method in the new syntax but a warning will be thrown. \code{distance = "logit"}, the old default, will still work in the new syntax, though \code{distance = "glm", link = "logit"} is preferred (note that these are the default settings and don't need to be made explicit).
}
\examples{
data("lalonde")
# Linearized probit regression PS:
m.out1 <- matchit(treat ~ age + educ + race + married +
                    nodegree + re74 + re75, data = lalonde,
                  distance = "glm", link = "linear.probit")
\dontshow{if (requireNamespace("CBPS", quietly = TRUE)) \{}
# GAM logistic PS with smoothing splines (s()):
m.out2 <- matchit(treat ~ s(age) + s(educ) + race + married +
                    nodegree + re74 + re75, data = lalonde,
                  distance = "gam")
summary(m.out2$model)
\dontshow{\}; if (requireNamespace("CBPS", quietly = TRUE)) \{}
# CBPS for ATC matching w/replacement, using the just-
# identified version of CBPS (setting method = "exact"):
m.out3 <- matchit(treat ~ age + educ + race + married +
                    nodegree + re74 + re75, data = lalonde,
                  distance = "cbps", estimand = "ATC",
                  distance.options = list(method = "exact"),
                  replace = TRUE)
\dontshow{\}}
# Mahalanobis distance matching - no PS estimated
m.out4 <- matchit(treat ~ age + educ + race + married +
                    nodegree + re74 + re75, data = lalonde,
                  distance = "mahalanobis")

m.out4$distance #NULL

# Mahalanobis distance matching with PS estimated
# for use in a caliper; matching done on mahvars
m.out5 <- matchit(treat ~ age + educ + race + married +
                    nodegree + re74 + re75, data = lalonde,
                  distance = "glm", caliper = .1,
                  mahvars = ~ age + educ + race + married +
                                nodegree + re74 + re75)

summary(m.out5)

# User-supplied propensity scores
p.score <- fitted(glm(treat ~ age + educ + race + married +
                        nodegree + re74 + re75, data = lalonde,
                      family = binomial))

m.out6 <- matchit(treat ~ age + educ + race + married +
                    nodegree + re74 + re75, data = lalonde,
                  distance = p.score)
}