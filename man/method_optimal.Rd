\name{method_optimal}
\alias{method_optimal}

\title{Optimal Pair Matching}

\description{
In \code{\link{matchit}}, setting \code{method = "optimal"} performs optimal pair matching. The matching is optimal in the sense that that sum of the absolute pairwise distances in the matched sample are as small as possible. The method relies on and is a wrapper for the \code{\link[optmatch]{pairmatch}} function in the \pkg{optmatch} package. Advantages of optimal pair matching include that the matching order is not required to be specified and it is less likely that extreme within-pair distances will be large, unlike with nearest neighbor matching. Generally, however, as a subset selection method, optimal pair matching tends to perform similarly to nearest neighbor matching in that similar subsets of units will be selected to be matched.

This page details the allowable arguments with \code{method = "optmatch"}. See \code{\link{matchit}} for an explanation of what each argument means in a general context and how it can be specified.

Below is how \code{matchit} is used for optimal pair matching:
\preformatted{
matchit(formula, data = NULL, method = "optimal",
        distance = "glm", link = "logit",
        distance.options = list(), estimand = "ATT",
        exact = NULL, mahvars = NULL, discard = "none",
        reestimate = FALSE, ratio = 1, verbose = FALSE,
        ...)
}
}
\arguments{
  \item{formula}{
a two-sided \code{\link{formula}} object containing the treatment and covariates to be used in creating the distance measure used in the matching. This formula will be supplied to the functions that estimate the distance measure.
}
  \item{data}{
a data frame containing the variables named in \code{formula}. If not found in \code{data}, the variables will be saught in the environment.
}
  \item{method}{
set here to \code{"optimal"}.
}
  \item{distance}{
the distance measure to be used. See \code{\link{distance}} for allowable options. When set to \code{"mahalanobis"}, optimal Mahalanobis distance matching will be performed on the variables named in \code{formula}.
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
for which variables exact matching should take place. Exact matching is processed using \code{\link[optmatch]{exactMatch}} in \pkg{optmatch}.
}
  \item{mahvars}{
for which variables Mahalanobis matching should take place when a distance measure other than \code{"mahalanobis"} is used (e.g., to discard units for common support). If specified, the distance measure will not be used in matching.
}
  \item{discard}{
a string containing a method for discarding units outside a region of common support. Only allowed when \code{distance} is not \code{"mahalanobis"}.
}
  \item{reestimate}{
if \code{discard} is not \code{"none"}, whether to re-estimate the propensity score in the remaining sample prior to matching.
}
  \item{ratio}{
how many control units should be matched to each treated unit for k:1 matching. Should be a single integer value.
}
  \item{verbose}{
\code{logical}; whether information about the matching process should be printed to the console. What is printed depends on the matching method. Default is \code{FALSE} for no printing other than warnings.
}
  \item{\dots}{
additional arguments passed to \code{\link[optmatch]{pairmatch}} in \pkg{optmatch}. Only the \code{tols} argument, which is eventually passed to \code{\link[optmatch]{fullmatch}}, make sense to specify.
}

The arguments \code{replace}, \code{caliper}, and \code{m.order} are ignored with a warning.
}
\note{
Calipers may eventually be compatible with this method, but due to a bug in \pkg{optmatch} (version 0.9-13), they are currently disabled.
}
\section{Outputs}{
All outputs described in \code{\link{matchit}} are returned with \code{method = "optimal"}.
}
\seealso{
\code{\link{matchit}} for a detailed explanation of the inputs and outputs of a call to \code{matchit}.

\code{\link[optmatch]{pairmatch}} in \pkg{optmatch}, which is the workhorse.

\code{\link{method_full}} for optimal full matching, of which optimal pair matching is a special case, and which relies on similar machinery.
}
\details{
Mahalanobis distance matching can be done one of two ways:

1) If a propensity score needs to be estimated for common support with \code{discard}, \code{distance} should be whatever method is used to estimate the propensity score or a vector of distance measures, i.e., it should not be \code{"mahalanobis"}. Use \code{mahvars} to specify the variables used to create the Mahalanobis distance. For example, to perform Mahalanobis after discarding units outside the common support of the propensity score in both groups, the following could be run:

\code{matchit(treat ~ X1 + X2 + X3, method = "optimal", distance = "glm",
              discard = "both", mahvars = ~ X1 + X2)}

With this code, \code{X1}, \code{X2}, and \code{X3} are used to estimate the propensity score (using the \code{"glm"} method, which by default is logistic regression), which is used to identify the common support. The actual matching occurs on the Mahalanobis distance computed only using \code{X1} and \code{X2}, which are supplied to \code{mahvars}. The estimated propensity scores will be included in the \code{distance} component of the \code{matchit} output.

2) If no propensity score needs to be estimated, \code{distance} should be set to \code{"mahalanobis"}, and Mahalanobis distance matching will occur on all the variables in \code{formula}. Arguments to \code{discard} and \code{mahvars} will be ignored. For example, to perform simple Mahalanobis distance matching, the following could be run:

\code{matchit(treat ~ X1 + X2, method = "optimal", distance = "mahalanobis")}

With this code, the Mahalanobis distance is computed using \code{X1} and \code{X2}, and matching occurs on this distance. The \code{distance} component of the \code{matchit} output will be empty.

\subsection{Estimand}{
The \code{estimand} argument controls whether control units are selected to be matched with treated units (\code{estimand = "ATT"}) or treated units are selected to be matched with control units (\code{estimand = "ATC"}). The "focal" group (e.g., the treated units for the ATT) is typically made to be the smaller treatmentgroup, and a warning will be thrown if it is not set that way unless \code{replace = TRUE}. Setting \code{estimand = "ATC"} is equivalent to swapping all treated and control labels for the treatment variable. When \code{estimand = "ATC"}, the \code{match.matrix} component of the output will have the names of the control units as the rownames and be filled with the names of the matched treated units (opposite to when \code{estimand = "ATT"}). Note that the argument supplied to \code{estimand} doesn't necessarily correspond to the estimand actually targeted; it is merely a switch to trigger which treatment group is considered "focal".
}

}
\references{
In a manuscript, be sure to cite the following paper if using \code{matchit} with \code{method = "optimal"}:

Hansen, B. B., & Klopfer, S. O. (2006). Optimal Full Matching and Related Designs via Network Flows. Journal of Computational and Graphical Statistics, 15(3), 609–627. \doi{10.1198/106186006X137047}

For example, a sentence might read:

\emph{Optimal pair matching was performed using the MatchIt package (Ho, Imai, King, & Stuart, 2011) in R, which calls functions from the optmatch package (Hansen & Klopfer, 2006).}
}
\examples{
data("lalonde")

# 1:1 optimal PS matching with exact matching on race
m.out1 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75, data = lalonde,
                  method = "optimal", exact = ~race)
m.out1
summary(m.out1)

#2:1 optimal Mahalanobis distance matching
m.out2 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75, data = lalonde,
                  method = "optimal", distance = "mahalanobis",
                  ratio = 2)
m.out2
summary(m.out2)
}