\name{distance}
\alias{distance}

\title{Propensity scores and other distance measures}

\description{
Several matching methods require or can involve the distance between treated and control units. Options include the Mahalanobis distance, propensity score distance, or distance between user-supplied values. Propensity scores are also used for common support via the \code{discard} options and for defining calipers. This page documents the options that can be supplied to the \code{distance} argument to \fun{matchit}.

There are four ways to specify the \code{distance} argument: 1) as the string \code{"mahalanobis"}, 2) as a string containing the name of a method for estimating propensity scores, 3) as a vector of values whose pairwise differences define the distance between units, or 4) as a distance matrix containing all pairwise differences.

When \code{distance} is specified as one of the allowed strings (described below) other than \code{"mahalanobis"}, a propensity score is estimated using the variables in \code{formula} and the method corresponding to the given argument. This propensity score can be used to compute the distance between units as the absolute difference between the propensity scores of pairs of units. In this respect, the propensity score is more like a "position" measure than a distance measure, since it is the pairwise differences that form the distance rather than the propensity scores themselves. Still, this naming convention is used to reflect their primary purpose without committing to the status of the estimated values as propensity scores, since transformations of the scores are allowed and user-supplied values that are not propensity scores can also be supplied (detailed below). Propensity scores can also be used to create calipers and common support restrictions, whether or not they are used in the actual distance measure used in the matching, if any.

In addition to the \code{distance} argument, two other arguments can be specified that relate to the estimation and manipulation of the propensity scores. The \code{link} argument allows for different links to be used in models that require them such as generalized linear models, for which the logit and probit links are allowed, among others. In addition to specifying the link, the \code{link} argument can be used to specify whether the propensity score or the linearized version of the propensity score should be used; by specifying \code{link = "linear.{link}"}, the linearized version will be used.

The \code{distance.options} argument can also be specified, which should be a list of values passed to the propensity score-estimating function, for example, to choose specific options or tuning parameters for the estimation method. If \code{formula}, \code{data}, or \code{verbose} are not supplied to \code{distance.options}, the corresponding arguments from \code{matchit()} will be automatically supplied. See the Examples for demonstrations of the uses of \code{link} and \code{distance.options}. When \code{s.weights} is supplied in the call to \code{matchit()}, it will automatically be passed to the propensity score-estimating function as the \code{weights} argument unless otherwise described below.
}

\section{Allowable options}{
Below are the allowed options for \code{distance}:

\describe{
  \item{\code{"glm"}}{
The propensity scores are estimated using a generalized linear model (e.g., logistic regression). The \code{formula} supplied to \code{matchit()} is passed directly to \fun{glm}, and \fun{predict.glm} is used to compute the propensity scores. The \code{link} argument can be specified as a link function supplied to \fun{binomial}, e.g., \code{"logit"}, which is the default. When \code{link} is prepended by \code{"linear."}, the linear predictor is used instead of the predicted probabilities. \code{distance = "glm"} with \code{link = "logit"} (logistic regression) is the default in \code{matchit()}.
}
  \item{\code{"gam"}}{
The propensity scores are estimated using a generalized additive model. The \code{formula} supplied to \code{matchit()} is passed directly to \pkgfun2{mgcv}{gam}{mgcv::gam}, and \pkgfun2{mgcv}{predict.gam}{mgcv::predict.gam} is used to compute the propensity scores. The \code{link} argument can be specified as a link function supplied to \fun{binomial}, e.g., \code{"logit"}, which is the default. When \code{link} is prepended by \code{"linear."}, the linear predictor is used instead of the predicted probabilities. Note that unless the smoothing functions \pkgfun{mgcv}{s}, \pkgfun{mgcv}{te}, \pkgfun{mgcv}{ti}, or \pkgfun{mgcv}{t2} are used in \code{formula}, a generalized additive model is identical to a generalized linear model and will estimate the same propensity scores as \code{glm}. See the documentation for \pkgfun2{mgcv}{gam}{mgcv::gam}, \pkgfun2{mgcv}{formula.gam}{mgcv::formula.gam}, and \pkgfun2{mgcv}{gam.models}{mgcv::gam.models} for more information on how to specify these models. Also note that the formula returned in the \code{matchit()} output object will be a simplified version of the supplied formula with smoothing terms removed (but all named variables present).
}
  \item{\code{"gbm"}}{
The propensity scores are estimated using a generalized boosted model. The \code{formula} supplied to \code{matchit()} is passed directly to \pkgfun2{gbm}{gbm}{gbm::gbm}, and \pkgfun2{gbm}{predict.gbm}{gbm::predict.gbm} is used to compute the propensity scores. The optimal tree is chosen using 5-fold cross-validation by default, and this can be changed by supplying an argument to \code{method} to \code{distance.options}; see \pkgfun2{gbm}{gbm.perf}{gbm::gbm.perf} for details. The \code{link} argument can be specified as \code{"linear"} to use the linear predictor instead of the predicted probabilities. No other links are allowed. The tuning parameter defaults differ from \code{gbm::gbm()}; they are as follows: \code{n.trees = 1e4}, \code{interaction.depth = 3}, \code{shrinkage = .01}, \code{bag.fraction = 1}, \code{cv.folds = 5}, \code{keep.data = FALSE}. These are the same defaults as used in \pkg{WeightIt} and \pkg{twang}, except for \code{cv.folds} and \code{keep.data}. Note this is not the same use of generalized boosted modeling as in \pkg{twang}; here, the number of trees is chosen based on cross-validation or out-of-bag error, rather than based on optimizing balance. \pkg{twang} should not be cited when using this method to estimate propensity scores.
}
  \item{\code{"lasso"}, \code{"ridge"}, \code{"elasticnet"}}{
The propensity scores are estimated using a lasso, ridge, or elastic net model, respectively. The \code{formula} supplied to \code{matchit()} is processed with \fun{model.matrix} and passed to \pkgfun2{glmnet}{cv.glmnet}{glmnet::cv.glmnet}, and \pkgfun2{glmnet}{predict.cv.glmnet}{glmnet::predict.cv.glmnet} is used to compute the propensity scores. The \code{link} argument can be specified as a link function supplied to \fun{binomial}, e.g., \code{"logit"}, which is the default. When \code{link} is prepended by \code{"linear."}, the linear predictor is used instead of the predicted probabilities. When \code{link = "log"}, a Poisson model is used. For \code{distance = "elasticnet"}, the \code{alpha} argument, which controls how to prioritize the lasso and ridge penalties in the elastic net, is set to .5 by default and can be changed by supplying an argument to \code{alpha} in \code{distance.options}. For \code{"lasso"} and \code{"ridge"}, \code{alpha} is set to 1 and 0, respectively, and cannot be changed. The \code{cv.glmnet()} defaults are used to select the tuning parameters and generate predictions and can be modified using \code{distance.options}. If the \code{s} argument is passed to \code{distance.options}, it will be passed to \code{predict.cv.glmnet()}. Note that because there is a random component to choosing the tuning parameter, results will vary across runs unless a \link[=set.seed]{seed} is set.
}
  \item{\code{"rpart"}}{
The propensity scores are estimated using a classification tree. The \code{formula} supplied to \code{matchit()} is passed directly to \pkgfun2{rpart}{rpart}{rpart::rpart}, and \pkgfun2{rpart}{predict.rpart}{rpart::predict.rpart} is used to compute the propensity scores. The \code{link} argument is ignored, and predicted probabilities are always returned as the distance measure.
}
  \item{\code{"randomforest"}}{
The propensity scores are estimated using a random forest. The \code{formula} supplied to \code{matchit()} is passed directly to \pkgfun2{randomForest}{randomForest}{randomForest::randomForest}, and \pkgfun2{randomForest}{predict.randomForest}{randomForest::predict.randomForest} is used to compute the propensity scores. The \code{link} argument is ignored, and predicted probabilities are always returned as the distance measure.% When \code{s.weights} is supplied to \code{matchit()}, it will not be passed to \code{randomForest} because \code{randomForest} does not accept weights.
}
  \item{\code{"nnet"}}{
The propensity scores are estimated using a single-hidden-layer neural network. The \code{formula} supplied to \code{matchit()} is passed directly to \pkgfun2{nnet}{nnet}{nnet::nnet}, and \fun{fitted} is used to compute the propensity scores. The \code{link} argument is ignored, and predicted probabilities are always returned as the distance measure. An argument to \code{size} must be supplied to \code{distance.options} when using \code{method = "nnet"}.
}
  \item{\code{"cbps"}}{
The propensity scores are estimated using the covariate balancing propensity score (CBPS) algorithm, which is a form of logistic regression where balance constraints are incorporated to a generalized method of moments estimation of of the model coefficients. The \code{formula} supplied to \code{matchit()} is passed directly to \pkgfun2{CBPS}{CBPS}{CBPS::CBPS}, and \code{\link{fitted}} is used to compute the propensity scores. The \code{link} argument can be specified as \code{"linear"} to use the linear predictor instead of the predicted probabilities. No other links are allowed. The \code{estimand} argument supplied to \code{matchit()} will be used to select the appropriate estimand for use in defining the balance constraints, so no argument needs to be supplied to \code{ATT} in \code{CBPS}.
}
  \item{\code{"bart"}}{
The propensity scores are estimated using Bayesian additive regression trees (BART). The \code{formula} supplied to \code{matchit()} is passed directly to \pkgfun2{dbarts}{bart}{dbarts::bart2}, and \pkgfun2{dbarts}{bart}{dbarts::fitted} is used to compute the propensity scores. The \code{link} argument can be specified as \code{"linear"} to use the linear predictor instead of the predicted probabilities. When \code{s.weights} is supplied to \code{matchit()}, it will not be passed to \code{bart2} because the \code{weights} argument in \code{bart2} does not correspond to sampling weights.
}
  \item{\code{"mahalanobis"}}{
No propensity scores are estimated. Rather than using the propensity score difference as the distance between units, the Mahalanobis distance is used instead. See \fun{mahalanobis} for details on how it is computed. The Mahalanobis distance is always computed using all the variables in \code{formula}. With this specification, calipers and common support restrictions cannot be used and the \code{distance} component of the output object will be empty because no propensity scores are estimated. The \code{link} and \code{distance.options} arguments are ignored. See individual methods pages for whether the Mahalanobis distance is allowed and how it is used. Sometimes this setting is just a placeholder to indicate that no propensity score is to be estimated (e.g., with \code{method = "genetic"}). To perform Mahalanobis distance matching \emph{and} estimate propensity scores to be used for a purpose other than matching, the \code{mahvars} argument should be used along with a different specification to \code{distance}. See the individual matching method pages for details on how to use \code{mahvars}.
}
}

\code{distance} can also be supplied as a numeric vector whose values will be taken to function like propensity scores; their pairwise difference will define the distance between units. This might be useful for supplying propensity scores computed outside \code{matchit()} or resupplying \code{matchit()} with propensity scores estimated before without having to recompute them. \code{distance} can also be supplied as a matrix whose values represent the pairwise distances between units. The matrix should either be a square, with a row and column for each unit (e.g., as the output of a call to \code{as.matrix(\link{dist}(.))}), or have as many rows as there are treated units and as many columns as there are control units (e.g., as the output of a call to \pkgfun2{optmatch}{match_on}{optmatch::match_on}). Distance values of \code{Inf} will disallow the corresponding units to be matched. When \code{distance} is a supplied as a numeric vector or matrix, \code{link} and \code{distance.options} are ignored.
}
\section{Outputs}{
When specifying an argument to \code{distance} that estimates a propensity score, the output of the function called to estimate the propensity score (e.g., the \code{glm} object when \code{distance = "glm"}) will be included in the \code{matchit()} output object in the \code{model} component. When \code{distance} is anything other than \code{"mahalanobis"} and not a matrix, the estimated or supplied distance measures will be included in the \code{matchit()} output object in the \code{distance} component.
}
\note{
In versions of \emph{MatchIt} prior to 4.0.0, \code{distance} was specified in a slightly different way. When specifying arguments using the old syntax, they will automatically be converted to the corresponding method in the new syntax but a warning will be thrown. \code{distance = "logit"}, the old default, will still work in the new syntax, though \code{distance = "glm", link = "logit"} is preferred (note that these are the default settings and don't need to be made explicit).
}
\examples{
data("lalonde")
# Linearized probit regression PS:
m.out1 <- matchit(treat ~ age + educ + race + married +
                    nodegree + re74 + re75, data = lalonde,
                  distance = "glm", link = "linear.probit")
\dontshow{if (requireNamespace("mgcv", quietly = TRUE)) \{}
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

# User-supplied distance matrix using optmatch::match_on()
\dontshow{if (requireNamespace("optmatch", quietly = TRUE)) \{}
dist_mat <- optmatch::match_on(
              treat ~ age + educ + race + nodegree +
                married + re74 + re75, data = lalonde,
              method = "rank_mahalanobis")

m.out7 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75, data = lalonde,
                  distance = dist_mat)
\dontshow{\}}}
