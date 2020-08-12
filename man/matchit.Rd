\name{matchit}
\alias{matchit}
\alias{MatchIt}
\alias{print.matchit}
\title{
Matching for Causal Inference
}
\description{
\code{matchit} is the main function of \pkg{MatchIt} and performs pairing, subset selection, and subclassification with the aim of creating treatment and control groups balanced on included covariates. MatchIt implements the suggestions of Ho, Imai, King, and Stuart (2007) for improving parametric statistical models by preprocessing data with nonparametric matching methods.  MatchIt implements a wide range of sophisticated matching methods, making it possible to greatly reduce the dependence of causal inferences on hard-to-justify, but commonly made, statistical modeling assumptions. The software also easily fits into existing research practices since, after preprocessing with MatchIt, researchers can use whatever parametric model they would have used without MatchIt, but produce inferences with substantially more robustness and less sensitivity to modeling assumptions.

This page documents the overall use of \code{matchit}, but for specifics of how \code{matchit} works with individual matching methods, see the individual pages linked in the Details section below.
}
\usage{
matchit(formula, data = NULL, method = "nearest",
        distance = "glm", link = "logit",
        distance.options = list(), estimand = "ATT",
        exact = NULL, mahvars = NULL, discard = "none",
        reestimate = FALSE, replace = FALSE,
        m.order = NULL, caliper = NULL, std.caliper = TRUE,
        ratio = 1, verbose = FALSE, ...)

\method{print}{matchit}(x, ...)

}

\arguments{
  \item{formula}{
a two-sided \code{\link{formula}} object containing the treatment and covariates to be used in creating the distance measure used in the matching. This formula will be supplied to the functions that estimate the distance measure. The formula should be specified as \code{A ~ X1 + X2 + ...} where \code{A} represents the treatment variable and \code{X1} and \code{X2} are covariates.
}
  \item{data}{
a data frame containing the variables named in \code{formula}. If not found in \code{data}, the variables will be saught in the environment.
}
  \item{method}{
the mathcing method to be used. The allowed methods are \code{\link[=method_nearest]{"nearest"}} for nearest neighbor matching (on the propensity score by default), \code{\link[=method_optimal]{"optimal"}} for optimal pair matching, \code{\link[=method_full]{"full"}} for optimal full matching, \code{\link[=method_genetic]{"genetic"}} for genetic matching, \code{\link[=method_cem]{"cem"}} for coarsened exact matching, \code{\link[=method_exact]{"exact"}} for exact matching, and \code{\link[=method_subclass]{"subclass"}} for subclassification. See the linked pages for each method for more details on what these methods do, how the arguments below are used by each on, and what additional arguments are allowed.
}
  \item{distance}{
the distance measure to be used. Can be either a string containing the name of a distance measure or a vector of already-computed distance measures. The distance measures should be values whose pairwise difference is the distance between two units, e.g., propensity scores for propensity score matching. See \code{\link{distance}} for allowable options. The default is \code{"glm"} for propensity scores estimated with logistic regression using \code{\link{glm}}. Ignored for some methods; see individual methods pages for how the distance measure is used.
}
  \item{link}{
when \code{distance} is specified as a string, an additional argument controlling the link function used in estimating the distance measure. Allowable options depend on the specific \code{distance} value specified. See \code{\link{distance}} for allowable options with each option. The default is \code{"logit"}, which, along with \code{distance = "glm"}, identifies the default measure as logistic regression propensity scores.
}
  \item{distance.options}{
a named list containing additional arguments supplied to the function that estimates the distance measure as determined by the argument to \code{distance}.
}
  \item{estimand}{
a string containing the name of the target estimand desired. Can be one of \code{"ATT"} or \code{"ATC"}. Some methods accept \code{"ATE"} as well. Default is \code{"ATT"}. See Details and the individual methods pages for information on how this argument is used.
}
  \item{exact}{
for methods that allow it, for which variables exact matching should take place when a matching method other than \code{"exact"} or \code{"cem"} is used. Can be specified as a string containing the names of variables in \code{data} to be used or a one-sided formula with the desired variables on the right-hand side (e.g., \code{~ X3 + X4}). See the individual methods pages for information on whether and how this argument is used.
}
  \item{mahvars}{
for methods that allow it, for which variables Mahalanobis matching should take place when a distance measure other than \code{"mahalanobis"} is used. Usually used to perform Mahalnobis distance matching within propensity score calipers, where the propensity scores are computed using \code{formula} and \code{distance}. Can be specified as a string containing the names of variables in \code{data} to be used or a one-sided formula with the desired variables on the right-hand side (e.g., \code{~ X3 + X4}). See the individual methods pages for information on whether and how this argument is used.
}
  \item{discard}{
a string containing a method for discarding units outside a region of common support. When a propensity score is estimated or supplied to \code{distance}, the options are \code{"none"}, \code{"treated"}, \code{"control"}, or \code{"both"}. For \code{"none"}, no units are discarded for common support. Otherwise, units whose propensity scores fall outside the corresponding region are discarded. Can also be a \code{logical} vector where \code{TRUE} indicates the unit is to be discarded. Default is \code{"none"} for no common support. See Details.
}
  \item{reestimate}{
if \code{discard} is not \code{"none"}, whether to re-estimate the propensity score in the remaining sample. Default is \code{FALSE} to use the propensity scores estimated in the original sample.
}
  \item{replace}{
for methods that allow it, whether matching should be done with replacement (\code{TRUE}), where control units are allowed to be matched to several treated units, or without replacement (\code{FALSE}), where controls units can only be matched to one treated unit each. See the individual methods pages for information on whether and how this argument is used. Default is \code{FALSE} for matching without replacement.
}
  \item{m.order}{
for methods that allow it, the order that the matching takes place. Allowable otpions depend on the matching method but include \code{"largest"}, where matching takes place in descending order of distance measures; \code{"smallest"}, where matching takes place in ascending order of distance measures; \code{"random"}, where matching takes place in a random order; and \code{"data"} where matching takes place based on the order of units in the data. When \code{m.order = "random"}, results may differ across different runs of the same code unless a seed is set and specified with \code{\link{set.seed}}. See the individual methods pages for information on whether and how this argument is used. The default of \code{NULL} corresponds to \code{"largest"} when a propensity score is estimated or supplied and \code{"data"} otherwise.
}
  \item{caliper}{
for methods that allow it, the width(s) of the caliper(s) to use in matching. Should be a numeric vector with each value named according to the variable to which the caliper applies. To apply to the distance measure, the value should be unnamed. See the individual methods pages for information on whether and how this argument is used. The default is \code{NULL} for no caliper.
}
  \item{std.caliper}{
\code{logical}; when a caliper is specified, whether the the caliper is in standard deviation units (\code{TRUE}) or raw units (\code{FALSE}). Can either be of length 1, applying to all calipers, or of length equal to the length of \code{caliper}. Default is \code{TRUE}.
}
  \item{ratio}{
for methods that allow it, how many control units should be matched to each treated unit in k:1 matching. Should be a single integer value. See the individual methods pages for information on whether and how this argument is used. The default is 1 for 1:1 matching.
}
  \item{verbose}{
\code{logical}; whether information about the matching process should be printed to the console. What is printed depends on the matching method. Default is \code{FALSE} for no printing other than warnings.
}
  \item{\dots}{
additional arguments passed to the functions used in the matching process. See the individual methods pages for information on what additional arguments are allowed for each method. Ignored for \code{print}.
}
  \item{x}{
a \code{matchit} object.
}
}
\details{
Details for the various matching methods can be found at the following help pages:
\itemize{
  \item{\code{\link{method_nearest}} for nearest neighbor matching}
  \item{\code{\link{method_optimal}} for optimal pair matching}
  \item{\code{\link{method_full}} for optimal full matching}
  \item{\code{\link{method_genetic}} for genetic matching}
  \item{\code{\link{method_cem}} for coarsened exact matching}
  \item{\code{\link{method_exact}} for exact matching}
  \item{\code{\link{method_subclass}} for subclassification}
}
The pages contain information on what the method does, which of the arguments above are allowed with them and how they are interpreted, and what additional arguments can be supplied to further tune the method. Note that the default method with no arguments supplied other than \code{formula} and \code{data} is 1:1 nearest neighbor matching without replacement on a propensity score estimated using a logistic regression of the treatment on the covariates. This is not the same default offered by other matching programs, such as those in \pkg{Matching}, \code{teffects} in Stata, or \code{PROC PSMATCH} in SAS, so care should be taken if trying to replicate the results of those programs.

See the \code{\link{distance}} argument for details on the several ways to specify the \code{distance} and \code{link} arguments to estimate propensity scores and create distance measures.

When the treatment variable is not a \code{0/1} variable, it will be coerced to one and returned as such in the \code{matchit} output (see section Value, below). The following rules are used: 1) if \code{0} is one of the values, it will be considered the control and the other value the treated; otherwise, 2) if the variable is a factor, \code{levels(treat)[1]} will be considered control and the other variable the treated; otherwise, 3) \code{sort(unique(treat))[1]} will be considered control and the other variable the treated. It is safest to ensure the treatment variable is a \code{0/1} variable.

The \code{discard} option implements a common support restriction. It can only be used when a distance measure is estimated or supplied, i.e., when \code{distance} is specified as something other than \code{"mahalanobis"}, and is ignored for some matching methods. When specified as \code{"treated"}, treated units whose distance measure is outside the range of distance measures of the control units will be discarded. When specified as \code{"control"}, control units whose distance measure is outside the range of distance measures of the treated units will be discarded. When specified as \code{"both"}, treated and controls units whose distance measure is outside the intersection of the range of distance measures of the treated units and the range of distance measures of the control units will be discarded. When \code{reestimate = TRUE} and \code{distance} corresponds to a propensity score-estimating function, the propensity scores are re-estimated in the remaining units prior to being used for matching or calipers.

Caution should be used when interpreting effects estimated with various values of \code{estimand}. Setting \code{estimand = "ATT"} doesn't necessarily mean the average treatment effect in the treated is being estimated; it just means that for matching methods, treated units will be untouched and given weights of 1 and control units will be matched to them (and the opposite for \code{estimand = "ATC"}). If a caliper is supplied or treated units are removed for common support or some other reason (e.g., lacking matches when using exact matching), the actual estimand targeted is not the ATT but the treatment effect in the matched sample. The argument to \code{estimand} simply triggers which units are matched to which, and for subclassification-based methods (exact matching, CEM, full matching, and subclassification), determines the formula used to compute the subclassification weights.

\subsection{Computing Weights}{
Matching weights are computed in one of two ways depending on whether matching was done with replacement or not.

For matching without replacement, each unit is assigned to a subclass, which represents the pair they are a part of (in the case of K:1 matching) or the subclass they belong to (in the case of exact matching, coarsened exact matching, full matching, or subclassification). The formula for computing the weights depends on the argument supplied to \code{estimand}. A new subclass "propensity score" (\code{p}) is computed as the proportion of units in each subclass that are in the treated group, and all units in that subclass are assigned that propensity score. Weights are then computed using the standard formulas for inverse probability weights: for the ATT, weights are 1 for the treated units and \code{p/(1-p)} for the control units; for the ATC, weights are \code{(1-p)/p} for the treated units and 1 for the control units; for the ATE, weights are \code{1/p} for the treated units and \code{1/(1-p)} for the control units.

For matching with replacement, units are not assigned to unique subclasses. For the ATT, each treated unit gets a weight of 1. Each control unit is weighted as the sum of the inverse of the number of control units matched to the same treated unit across its matches. For example, if a control unit was matched to a treated unit that had two other control units matched to it, and that same control was matched to a treated unit that had one other control unit matched to it, the control unit in question would get a weight of 1/3 + 1/2 = 5/6. For the ATC, the same is true with the treated and control labels switched. The weights are computed using the \code{match.matrix} component of the \code{matchit} output object.

In each treatment group, weights are divided by the sum of the weights in that treatment group to make the weights sum to the number of units in that treatment group.
}

}
\value{
When \code{method} is something other than \code{"subclass"}, a \code{matchit} object with the following components:
\itemize{
  \item{\code{match.matrix} - a matrix containing the matches. The rownames correspond to the treated units and the values in each row are the names (or indices) of the control units matched to each treated unit. When treated units are matched to different numbers of control units (e.g., with exact matching or matching with a caliper), empty spaces will be filled with \code{NA}.}
  \item{\code{subclass} - an integer vector of matching pair/stratum membership for each unit. Unmatched units will have a value of \code{NA}. Not included when \code{replace = TRUE}.}
  \item{\code{weights} - a numeric vector of estimated matching weights. Unmatched units will have a weight of zero.}
  \item{\code{model} - the fit object of the model used to estimate propensity scores when \code{distance} is specified and not \code{"mahalanobis"} or a numeric vector. When \code{reestimate = TRUE}, this is the model estimated after discarding units.}
  \item{\code{X} - a data frame of covariates mentioned in \code{formula}, \code{exact}, and \code{mahvars}.}
  \item{\code{call} - the \code{matchit()} call.}
  \item{\code{info} - information on the matching method and distance measures used.}
  \item{\code{estimand} - the argument supplied to \code{estimand}.}
  \item{\code{formula} - the \code{formula} supplied.}
  \item{\code{treat} - a vector of treatment status converted to zeros (0) and ones (1) if not already in that format.}
  \item{\code{distance} - a vector of distance values (i.e., propensity scores) when \code{distance} is specified and not \code{"mahalanobis"}.}
  \item{\code{discarded} - a logical vector denoting whther each observation was discarded (\code{TRUE}) or not (\code{FALSE}) by the argument to \code{discard}.}
  \item{\code{exact} - a one-sided formula containing the variables, if any, supplied to \code{exact}.}
  \item{\code{mahvars} - a one-sided formula containing the variables, if any, supplied to \code{mahvars}.}
  \item{\code{nn} - a matrix of the sample sizes of the treated and control groups before and after matching.}
}

When \code{method = "subclass"}, a \code{matchit.subclass} object with the same components as above except that \code{match.matrix} is excluded and an additional component, \code{q.cut}, is included, containing a vector of the distance measure cutpoints used to define the subclasses. See \code{\link{method_subclass}} for details.
}
\references{
Ho, D. E., Imai, K., King, G., & Stuart, E. A. (2007). Matching as Nonparametric Preprocessing for Reducing Model Dependence in Parametric Causal Inference. Political Analysis, 15(3), 199–236. \doi{10.1093/pan/mpl013}
}
\author{
Daniel Ho \email{dho@law.stanford.edu}; Kosuke Imai \email{imai@harvard.edu}; Gary King \email{king@harvard.edu}; Elizabeth Stuart \email{estuart@jhsph.edu}

2020 update by Noah Greifer \email{noah.greifer@gmail.com}
}

\seealso{
\code{\link{summary.matchit}} for balance assessment after matching. \code{\link{plot.matchit}} for plots of covariate balance and propensity score overlap after matching.
}
\examples{
data("lalonde")

# Default: 1:1 NN PS matching w/o replacement

m.out1 <- matchit(treat ~ age + educ + race + nodegree +
                   married + re74 + re75, data = lalonde)
m.out1
summary(m.out1)

# 1:1 NN Mahalanobis distance matching w/ replacement and
# exact matching on married and race

m.out2 <- matchit(treat ~ age + educ + race + nodegree +
                   married + re74 + re75, data = lalonde,
                   distance = "mahalanobis", replace = TRUE,
                   exact = ~ married + race)
m.out2
summary(m.out2)

# 2:1 NN Mahalanobis distance matching within caliper defined
# by a probit pregression PS

m.out3 <- matchit(treat ~ age + educ + race + nodegree +
                   married + re74 + re75, data = lalonde,
                   distance = "glm", link = "probit",
                   mahvars = ~ age + educ + re74 + re75,
                   caliper = .1, ratio = 2)
m.out3
summary(m.out3)

# Optimal full PS matching within calipers on PS, age, and educ

m.out4 <- matchit(treat ~ age + educ + race + nodegree +
                   married + re74 + re75, data = lalonde,
                   method = "full",
                   caliper = c(.1, age = 2, educ = 1),
                   std.caliper = c(TRUE, FALSE, FALSE)
m.out4
summary(m.out4)

# Subclassification on a CBPS PS with 10 subclasses after
# discarding controls outside common support of PS

s.out1 <- matchit(treat ~ age + educ + race + nodegree +
                   married + re74 + re75, data = lalonde,
                   method = "subclass", distance = "cbps",
                   discard = "control", subclass = 10)
s.out1
summary(s.out1)

}