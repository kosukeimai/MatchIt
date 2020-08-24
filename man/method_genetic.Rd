\name{method_genetic}
\alias{method_genetic}
\title{
Genetic Matching
}
\description{
In \code{\link{matchit}}, setting \code{method = "genetic"} performs genetic matching. Genetic matching is a form of nearest neighbor matching where distances are computed as the generalized Mahalanobis distance, which is a generalization of the Mahalanobis distance with a scaling factor for each covariate that represents the importance of that covariate. A genetic algorithm is used to select the scaling factors. The scaling factors are chosen as those which maximize a criterion related to covariate balance, which can be chosen, but which by default is the smallest p-value in covariate balance tests among the covariates. This method relies on and is a wrapper for the \code{\link[Matching]{GenMatch}} and \code{\link[Matching]{Match}} functions in the \pkg{Matching} package, which uses \code{\link[rgenoud]{genoud}} from the \pkg{rgenoud} package to perform the optimization using the genetic algorithm.

This page details the allowable arguments with \code{method = "genmatch"}. See \code{\link{matchit}} for an explanation of what each argument means in a general context and how it can be specified.

Below is how \code{matchit} is used for genetic matching:
\preformatted{
matchit(formula, data = NULL, method = "genetic",
        distance = "glm", link = "logit",
        distance.options = list(), estimand = "ATT",
        exact = NULL, mahvars = NULL, discard = "none",
        reestimate = FALSE, replace = FALSE,
        m.order = NULL, caliper = NULL, ratio = 1,
        verbose = FALSE, ...)
}
}
\arguments{
   \item{formula}{
a two-sided \code{\link{formula}} object containing the treatment and covariates to be used in creating the distance measure used in the matching. This formula will be supplied to the functions that estimate the distance measure and is used to determine the covariates whose balance is to be optimized.
}
  \item{data}{
a data frame containing the variables named in \code{formula}. If not found in \code{data}, the variables will be saught in the environment.
}
  \item{method}{
set here to \code{"genetic"}.
}
  \item{distance}{
the distance measure to be used. See \code{\link{distance}} for allowable options. When set to \code{"mahalanobis"}, only the covariates in \code{formula} are supplied to the generalized Mahalanobis distance matrix to have their scaling factors chosen. Otherwise, the distance measure is included with the covariates in \code{formula} to be supplied to the generalized Mahalanobis distance matrix.
}
  \item{link}{
when \code{distance} is specified as a string and not \code{"mahalanobis"}, an additional argument controlling the link function used in estimating the distance measure. See \code{\link{distance}} for allowable options with each option.
}
  \item{distance.options}{
a named list containing additional arguments supplied to the function that estimates the distance measure as determined by the argument to \code{distance}.
}
  \item{estimand}{
a string containing the desired estimand. Allowable options include \code{"ATT"} and \code{"ATC"}. See Details.
}
  \item{exact}{
for which variables exact matching should take place.
}
  \item{mahvars}{
when a distance measure other than \code{"mahalanobis"} is used (e.g., for caliper matching or to discard units for common support), which covariates should be supplied to the generalized Mahalanobis distance matrix. If unspecified, all variables in \code{formula} will be supplied to the distance matrix. Use \code{mahvars} to only supply a subset. Even if \code{mahvars} is specified, balance will be optimized on all covariates in \code{formula}.
}
  \item{discard}{
a string containing a method for discarding units outside a region of common support. Only allowed when \code{distance} is not \code{"mahalanobis"}.
}
  \item{reestimate}{
if \code{discard} is not \code{"none"}, whether to re-estimate the propensity score in the remaining sample prior to matching.
}
  \item{replace}{
whether matching should be done with replacement.
}
  \item{m.order}{
the order that the matching takes place. The default for \code{distance = "mahalanobis"} is \code{"data"}. Otherwise, the default is \code{"largest"}. See \code{\link{matchit}} for allowable options.
}
  \item{caliper}{
the width(s) of the caliper(s) used for caliper matching. See Details and Examples.
}
  \item{std.caliper}{
\code{logical}; when calipers are specified, whether they are in standard deviation units (\code{TRUE}) or raw units (\code{FALSE}).
}
  \item{ratio}{
how many control units should be matched to each treated unit for k:1 matching. Should be a single integer value.
}
  \item{verbose}{
\code{logical}; whether information about the matching process should be printed to the console. When \code{TRUE}, output from \code{GenMatch} with \code{print.level = 2} will be displayed. Default is \code{FALSE} for no printing other than warnings.
}
  \item{\dots}{
additional arguments passed to \code{\link[Matching]{GenMatch}} in \pkg{Matching}. Potentially useful options include \code{pop.size}, \code{max.generations}, and \code{fit.func}. If \code{pop.size} is not specified, a warning from \pkg{Matching} will be thrown reminding you to change it. Note that the \code{ties} and \code{CommonSupport} arguments are set to \code{FALSE} and cannot be changed.
}
}
\section{Outputs}{
All outputs described in \code{\link{matchit}} are returned with \code{method = "genetic"}. When \code{repalce = TRUE}, the \code{subclass} component is omitted.
}
\details{
In genetic matching, covariates play three roles: 1) as the variables on which balance is optimized, 2) as the variables in the generalized Mahalanobis distance between units, and 3) in estimating the propensity score. Variables supplied to \code{formula} are always used for role (1), as the variables on which balance is optimized. When \code{distance} is not \code{"mahalanobis"}, the covariates are also used to estimate the propensity score (unless it is supplied). When \code{mahvars} is specified, the named variables will form the covariates that go into the distance matrix. Otherwise, the variables in \code{formula} along with the propensity score will go into the distance matrix. This leads to three ways to use \code{distance} and \code{mahvars} to perform the matching:

1) When \code{distance = "mahalanobis"}, no propensity score is estimated, and the covariates in \code{formula} are used to form the generalized Mahalanobis distance matrix. In this sense, \code{"mahalanobis"} signals that no propensity score is to be estimated and that the matching variables are those in \code{formula}, consistent with setting \code{distance = "mahalanobis"} with other methods.

2) When \code{distance} is not \code{"mahalanobis"} and \code{mahvars} \emph{is not} specified, the covariates in \code{formula} along with the propensity score are used to form the generalized Mahalanobis distance matrix. This is the default and most typical use of \code{method = "genetic"} in \pkg{MatchIt}.

3) When \code{distance} is not \code{"mahalanobis"} and \code{mahvars} \emph{is} specified, the covariates in \code{mahvars} are used to form the generalized Mahalanobis distance matrix. The covariates in \code{formula} are used to estimate the propensity score and have their balance optimized by the genetic algorithm. The propensity score is not included in the generalized Mahalanobis distance matrix.

When a caliper is specified, any variables mentioned in \code{caliper}, possibly including the distance measure, will be added to the matching variables used to form the generalized Mahalanobis distance matrix. This is because \pkg{Matching} doesn't allow for the separation of caliper variables and matching variables in genetic matching.

\subsection{Estimand}{
The \code{estimand} argument controls whether control units are selected to be matched with treated units (\code{estimand = "ATT"}) or treated units are selected to be matched with control units (\code{estimand = "ATC"}). The "focal" group (e.g., the treated units for the ATT) is typically made to be the smaller treatmentgroup, and a warning will be thrown if it is not set that way unless \code{replace = TRUE}. Setting \code{estimand = "ATC"} is equivalent to swapping all treated and control labels for the treatment variable. When \code{estimand = "ATC"}, the default \code{m.order} is \code{"smallest"}, and the \code{match.matrix} component of the output will have the names of the control units as the rownames and be filled with the names of the matched treated units (opposite to when \code{estimand = "ATT"}). Note that the argument supplied to \code{estimand} doesn't necessarily correspond to the estimand actually targeted; it is merely a switch to trigger which treatment group is considered "focal". Note that while \code{GenMatch()} and \code{Matching()} support the ATE as an estimand, \code{matchit()} only supports the ATT and ATC for genetic matching.
}
}
\references{
In a manuscript, be sure to cite the following papers if using \code{matchit} with \code{method = "genetic"}:

Diamond, A., & Sekhon, J. S. (2013). Genetic matching for estimating causal effects: A general multivariate matching method for achieving balance in observational studies. Review of Economics and Statistics, 95(3), 932–945. \doi{10.1162/REST_a_00318}

Sekhon, J. S. (2011). Multivariate and Propensity Score Matching Software with Automated Balance Optimization: The Matching package for R. Journal of Statistical Software, 42(1), 1–52. \doi{10.18637/jss.v042.i07}

For example, a sentence might read:

\emph{Genetic matching was performed using the MatchIt package (Ho, Imai, King, & Stuart, 2011) in R, which calls functions from the Matching package (Diamond & Sekhon, 2013; Sekhon, 2011).}
}

\seealso{
\code{\link{matchit}} for a detailed explanation of the inputs and outputs of a call to \code{matchit}.

\code{\link[Matching]{GenMatch}} and \code{\link[Matching]{Match}} in \pkg{Matching}, which do the work.
}
\examples{
data("lalonde")

# 1:1 genetic matching with PS as a covariate
m.out1 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75, data = lalonde,
                  method = "genetic")
m.out1
summary(m.out1)

# 2:1 genetic matching with replacement without PS
m.out2 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75, data = lalonde,
                  method = "genetic", replace = TRUE,
                  ratio = 2, distance = "mahalanobis")
m.out2
summary(m.out2)

# 1:1 genetic matching on just age, educ, re74, and re75
# within calipers on PS and educ (PS is not included in
# the matching variables)
m.out3 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75, data = lalonde,
                  method = "genetic",
                  mahvars = ~ age + educ + re74 + re75,
                  caliper = c(.05, educ = 2),
                  std.caliper = c(TRUE, FALSE))
m.out3
summary(m.out3)
}