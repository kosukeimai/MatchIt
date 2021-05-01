\name{method_nearest}
\alias{method_nearest}
\title{
Nearest Neighbor Matching
}
\description{
In \fun{matchit}, setting \code{method = "nearest"} performs greedy nearest neighbor matching. A distance is computed between each treated unit and each control unit, and, one by one, each treated unit is assigned a control unit as a match. The matching is "greedy" in the sense that there is no action taken to optimize an overall criterion; each match is selected without considering the other matches that may occur subsequently.

This page details the allowable arguments with \code{method = "nearest"}. See \fun{matchit} for an explanation of what each argument means in a general context and how it can be specified.

Below is how \code{matchit()} is used for nearest neighbor matching:
\preformatted{
matchit(formula,
        data = NULL,
        method = "nearest",
        distance = "glm",
        link = "logit",
        distance.options = list(),
        estimand = "ATT",
        exact = NULL,
        mahvars = NULL,
        antiexact = NULL,
        discard = "none",
        reestimate = FALSE,
        s.weights = NULL,
        replace = TRUE,
        m.order = NULL,
        caliper = NULL,
        ratio = 1,
        min.controls = NULL,
        max.controls = NULL,
        verbose = FALSE, ...)
}
}
\arguments{
   \item{formula}{
a two-sided \fun{formula} object containing the treatment and covariates to be used in creating the distance measure used in the matching.
}
  \item{data}{
a data frame containing the variables named in \code{formula}. If not found in \code{data}, the variables will be sought in the environment.
}
  \item{method}{
set here to \code{"nearest"}.
}
  \item{distance}{
the distance measure to be used. See \code{\link{distance}} for allowable options. Can be supplied as a distance matrix.
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
for which variables Mahalanobis distance matching should take place when a distance measure other than \code{"mahalanobis"} is used (e.g., for caliper matching or to discard units for common support). If specified, the distance measure will not be used in matching.
}
  \item{antiexact}{
for which variables ant-exact matching should take place.
}
  \item{discard}{
a string containing a method for discarding units outside a region of common support. Only allowed when \code{distance} is not \code{"mahalanobis"} and not a matrix.
}
  \item{reestimate}{
if \code{discard} is not \code{"none"}, whether to re-estimate the propensity score in the remaining sample prior to matching.
}
  \item{s.weights}{
the variable containing sampling weights to be incorporated into propensity score models and balance statistics.
}
  \item{replace}{
whether matching should be done with replacement.
}
  \item{m.order}{
the order that the matching takes place. The default for \code{distance = "mahalanobis"} or when \code{distance} is supplied as a matrix is \code{"data"}. Otherwise, the default is \code{"largest"}. See \fun{matchit} for allowable options.
}
  \item{caliper}{
the width(s) of the caliper(s) used for caliper matching. See Details and Examples.
}
  \item{std.caliper}{
\code{logical}; when calipers are specified, whether they are in standard deviation units (\code{TRUE}) or raw units (\code{FALSE}).
}
  \item{ratio}{
how many control units should be matched to each treated unit for k:1 matching. For variable ratio matching, see section "Variable Ratio Matching" in Details below.
}
  \item{min.controls, max.controls}{
for variable ratio matching, the minimum and maximum number of controls units to be matched to each treated unit. See section "Variable Ratio Matching" in Details below.
}
  \item{verbose}{
\code{logical}; whether information about the matching process should be printed to the console. When \code{TRUE}, a progress bar implemented using \emph{RcppProgress} will be displayed.
}
  \item{\dots}{
ignored.
}
}
\section{Outputs}{
All outputs described in \fun{matchit} are returned with \code{method = "nearest"}. When \code{replace = TRUE}, the \code{subclass} component is omitted.
}
\details{
\subsection{Mahalanobis Distance Matching}{
Mahalanobis distance matching can be done one of two ways:

1) If no propensity score needs to be estimated, \code{distance} should be set to \code{"mahalanobis"}, and Mahalanobis distance matching will occur on all the variables in \code{formula}. Arguments to \code{discard} and \code{mahvars} will be ignored, and a caliper can only be placed on named variables. For example, to perform simple Mahalanobis distance matching, the following could be run:

\preformatted{
matchit(treat ~ X1 + X2, method = "nearest",
        distance = "mahalanobis")
}

With this code, the Mahalanobis distance is computed using \code{X1} and \code{X2}, and matching occurs on this distance. The \code{distance} component of the \code{matchit()} output will be empty.

2) If a propensity score needs to be estimated for any reason, e.g., for common support with \code{discard} or for creating a caliper, \code{distance} should be whatever method is used to estimate the propensity score or a vector of distance measures, i.e., it should not be \code{"mahalanobis"}. Use \code{mahvars} to specify the variables used to create the Mahalanobis distance. For example, to perform Mahalanobis within a propensity score caliper, the following could be run:

\preformatted{
matchit(treat ~ X1 + X2 + X3, method = "nearest",
        distance = "glm", caliper = .25,
        mahvars = ~ X1 + X2)
}

With this code, \code{X1}, \code{X2}, and \code{X3} are used to estimate the propensity score (using the \code{"glm"} method, which by default is logistic regression), which is used to create a matching caliper. The actual matching occurs on the Mahalanobis distance computed only using \code{X1} and \code{X2}, which are supplied to \code{mahvars}. Units whose propensity score difference is larger than the caliper will not be paired, and some treated units may therefore not receive a match. The estimated propensity scores will be included in the \code{distance} component of the \code{matchit()} output. See Examples.

When sampling weights are supplied through the \code{s.weights} argument, the covariance matrix of the covariates used in the Mahalanobis distance is weighted by the sampling weights.

Mahalanobis distance matching can also be done by supplying a Mahalanobis distance matrix (e.g., the output of a call to \pkgfun2{optmatch}{match_on}{optmatch::match_on}) to the \code{distance} argument. This makes it straightforward to use the robust rank-based Mahalanobis distance available in \pkg{optmatch}.
}
\subsection{Estimand}{
The \code{estimand} argument controls whether control units are selected to be matched with treated units (\code{estimand = "ATT"}) or treated units are selected to be matched with control units (\code{estimand = "ATC"}). The "focal" group (e.g., the treated units for the ATT) is typically made to be the smaller treatment group, and a warning will be thrown if it is not set that way unless \code{replace = TRUE}. Setting \code{estimand = "ATC"} is equivalent to swapping all treated and control labels for the treatment variable. When \code{estimand = "ATC"}, the default \code{m.order} is \code{"smallest"}, and the \code{match.matrix} component of the output will have the names of the control units as the rownames and be filled with the names of the matched treated units (opposite to when \code{estimand = "ATT"}). Note that the argument supplied to \code{estimand} doesn't necessarily correspond to the estimand actually targeted; it is merely a switch to trigger which treatment group is considered "focal".
}

\subsection{Variable Ratio Matching}{
\code{matchit()} can perform variable ratio "extremal" matching as described by Ming and Rosenbaum (2000). This method tends to result in better balance than fixed ratio matching at the expense of some precision. When \code{ratio > 1}, rather than requiring all treated units to receive \code{ratio} matches, each treated unit is assigned a value that corresponds to the number of control units they will be matched to. These values are controlled by the arguments \code{min.controls} and \code{max.controls}, which correspond to \ifelse{latex}{\out{$\alpha$}}{\ifelse{html}{\out{&alpha;}}{alpha}} and \ifelse{latex}{\out{$\beta$}}{\ifelse{html}{\out{&beta;}}{beta}}, respectively, in Ming and Rosenbaum (2000), and trigger variable ratio matching to occur. Some treated units will receive \code{min.controls} matches and others will receive \code{max.controls} matches (and one unit may have an intermediate number of matches); how many units are assigned each number of matches is determined by the algorithm described in Ming and Rosenbaum (2000, p119). \code{ratio} controls how many total control units will be matched: \code{n1 * ratio} control units will be matched, where \code{n1} is the number of treated units, yielding the same total number of matched controls as fixed ratio matching does.

Variable ratio matching cannot be used with Mahalanobis distance matching or when \code{distance} is supplied as a matrix. The calculations of the numbers of control units each treated unit will be matched to occurs without consideration of \code{caliper} or \code{discard}. \code{ratio} does not have to be an integer but must be greater than 1 and less than \code{n0/n1}, where \code{n0} and \code{n1} are the number of control and treated units, respectively. Setting \code{ratio = n0/n1} performs a crude form of full matching where all control units are matched. If \code{min.controls} is not specified, it is set to 1 by default. \code{min.controls} must be less than \code{ratio}, and \code{max.controls} must be greater than \code{ratio}. See Examples below for an example of their use.
}

}
\references{
In a manuscript, you don't need to cite another package when using \code{method = "nearest"} because the matching is performed completely within \emph{MatchIt}. For example, a sentence might read:

\emph{Nearest neighbor matching was performed using the MatchIt package (Ho, Imai, King, & Stuart, 2011) in R.}
}

\seealso{
\fun{matchit} for a detailed explanation of the inputs and outputs of a call to \code{matchit()}.

\code{\link{method_optimal}} for optimal pair matching, which is similar to nearest neighbor matching except that an overall distance criterion is minimized.
}
\examples{
data("lalonde")

# 1:1 greedy NN matching on the PS
m.out1 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75, data = lalonde,
                  method = "nearest")
m.out1
summary(m.out1)

# 3:1 NN Mahalanobis distance matching with
# replacement within a PS caliper
m.out2 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75, data = lalonde,
                  method = "nearest", replace = TRUE,
                  mahvars = ~ age + educ + re74 + re75,
                  ratio = 3, caliper = .02)
m.out2
summary(m.out2)

# 1:1 NN Mahalanobis distance matching within calipers
# on re74 and re75 and exact matching on married and race
m.out3 <- matchit(treat ~ age + educ + re74 + re75, data = lalonde,
                  method = "nearest", distance = "mahalanobis",
                  exact = ~ married + race,
                  caliper = c(re74 = .2, re75 = .15))
m.out3
summary(m.out3)

# 2:1 variable ratio NN matching on the PS
m.out4 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75, data = lalonde,
                  method = "nearest", ratio = 2,
                  min.controls = 1, max.controls = 12)
m.out4
summary(m.out4)

# Some units received 1 match and some received 12
table(table(m.out4$subclass[m.out4$treat == 0]))
}
