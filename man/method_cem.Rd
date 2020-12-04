\name{method_cem}
\alias{method_cem}

\title{Coarsened Exact Matching}

\description{
In \code{\link{matchit}}, setting \code{method = "cem"} performs coarsened exact matching. With coarsened exact matching, covariates are coarsened into bins, and a complete cross of the coarsened covariates is used to form subclasses defined by each combination of the coarsened covariate levels. Any subclass that doesn't contain both treated and control units is discarded, leaving only subclasses containing treatment and control units that are exactly equal on the coarsened covariates. The coarsening process can be controlled by an algorithm or by manually specifying cutpoints and groupings. The benefits of coarsened exact matching are that the tradeoff between exact matching and approximate balancing can be managed to prevent discarding too many units, which can otherwise occur with exact matching.

This page details the allowable arguments with \code{method = "cem"}. See \code{\link{matchit}} for an explanation of what each argument means in a general context and how it can be specified.

Below is how \code{matchit} is used for coarsened exact matching:
\preformatted{
matchit(formula, data = NULL, method = "cem",
        estimand = "ATT", s.weights = NULL,
        verbose = FALSE, ...)
}
}
\arguments{
   \item{formula}{
a two-sided \code{\link{formula}} object containing the treatment and covariates to be used in creating the subclasses defined by a full cross of the coarsened covariate levels.
}
  \item{data}{
a data frame containing the variables named in \code{formula}. If not found in \code{data}, the variables will be sought in the environment.
}
  \item{method}{
set here to \code{"cem"}.
}
  \item{estimand}{
a string containing the desired estimand. Allowable options include \code{"ATT"}, \code{"ATC"}, and \code{"ATE"}. The estimand controls how the weights are computed; see the Computing Weights section at \code{\link{matchit}} for details. When \code{k2k = TRUE} (see below), \code{estimand} also controls how the matching is done.
}
  \item{s.weights}{
the variable containing sampling weights to be incorporated into balance statistics. These weights do not affect the matching process.
}
  \item{verbose}{
\code{logical}; whether information about the matching process should be printed to the console.
}
  \item{\dots}{
  additional arguments to control the matching process.
  \describe{
     \item{\code{cutpoints}}{
        a named list with an (optional) entry for each numeric variable to be matched on. Each element describes a way of coarsening the corresponding variable. They can be a vector of cutpoints that demarcate bins, a single number giving the number of bins, or a string corresponding to a method of computing the number of bins. Allowable strings include \code{"sturges"}, \code{"scott"}, and \code{"fd"}, which use the functions \pkgfun2{grDevices}{nclass}{nclass.Sturges}, \pkgfun2{grDevices}{nclass}{nclass.scott}, and \pkgfun2{grDevices}{nclass}{nclass.FD}, respectively. The default is \code{"sturges"} for variables that are not listed or if no argument is supplied.
     }
     \item{\code{grouping}}{
     a named list with an (optional) entry for each categorical variable to be matched on. Each element should itself be a list, and each entry of the sublist should be a vector containing levels of the variable that should be combined to form a single level. For example, if a variable \code{amount} had levels \code{"none"}, \code{"some"}, and \code{"a lot"}, one could enter \code{grouping = list(amount = list("none", c("some", "a lot")))}, which would group \code{"some"} and \code{"a lot"} into a single category. Any categorical variables not included in \code{grouping} will remain as they are in the data, which means exact matching, with no coarsening, will take place on these variables.
     }
     \item{\code{k2k}}{
     code{logical}; whether 1:1 matching should occur within the matched strata. If \code{TRUE} nearest neighbor matching without replacement will take place within each stratum, and any unmatched units will be dropped (e.g., if there are more treated than control units in the stratum, the treated units without a match will be dropped). The \code{k2k.method} argument controls how the distance between units is calculated.
     }
     \item{\code{k2k.method}}{
     \code{character}; how the distance between units should be calculated if \code{k2k = TRUE}. Allowable arguments include \code{NULL} (for random matching), \code{"mahalanobis"} (for Mahalanobis distance matching), or any allowable argument to \code{method} in \code{\link{dist}}. Matching will take place on scaled versions of the original (non-coarsened) variables. The default is \code{"mahalanobis"}.
     }
     \item{\code{mpower}}{
     if \code{k2k.method = "minkowski"}, the power used in creating the distance. This is passed to the \coe{p} argument of \code{\link{dist}}.
     }
  }
}
The arguments \code{distance} (and related arguments), \code{exact}, \code{mahvars}, \code{discard} (and related arguments), \code{replace}, \code{m.order}, \code{caliper} (and related arguments), and \code{ratio} are ignored with a warning.
}
\section{Outputs}{
All outputs described in \code{\link{matchit}} are returned with \code{method = "cem"} except for \code{match.matrix}. When \code{k2k = TRUE}, a \code{match.matrix} component with the matched pairs is also included.
}
\details{
If the coarsening is such that there are no exact matches with the coarsened variables, the \code{cutpoints} and \code{grouping} arguments can be used to modify the matching specification. Reducing the number of cutpoints or grouping some variable values together can make it easier to find matches. See Examples below. Removing variables can also help (but they will likely not be balanced unless highly correlated with the included variables). To take advantage of coarsened exact matching without failing to find any matches, the covariates can be manually coarsened outside of \code{matchit()} and then supplied to the \code{exact} argument in a call to \code{matchit()} with another matching method.

This method does not rely on the \pkg{cem} package, instead using code written for \pkg{MatchIt}, but its design is based on the original \pkg{cem} functions. Versions of \pkg{MatchIt} prior to 4.1.0 did rely on \pkg{cem}, so results may differ between versions. There are a few differences between the ways \pkg{MatchIt} and \pkg{cem} (and older versions of \pkg{MatchIt}) differ in executing coarsened exact matching, described below.
\itemize{
\item{In \pkg{MatchIt}, when a single number is supplied to \code{cutpoints}, it describes the number of bins; in \pkg{cem} it describes the number of cutpoints separating bins. The \pkg{MatchIt} method is closer to how \code{\link{hist}} processes breaks points.}
\item{When \code{cutpoints} are used, \code{"ss"} (for Shimazaki-Shinomoto's rule) can be used in \pkg{cem} but not in \pkg{MatchIt}.}
\item{When \code{k2k = TRUE}, \pkg{MatchIt} matches on the original variables (scaled), whereas \pkg{cem} matches on the coarsened variables. Because the variables are already exactly matched on the coarsened variables, matching in \code{cem} is equivalent to random matching within strata.}
\item{When \code{k2k = TRUE}, in \pkg{MatchIt} matched units are identified by pair membership, and the original stratum membership prior to 1:1 matching is discarded. In \pkg[cem}, pairs are not identified beyond the stratum the members are part of.}
\item{When \code{k2k = TRUE}, \code{k2k.method = "mahalanobis"} can be requested in \pkg{MatchIt} but not in \pkg{cem}.}
}

Setting \code{k2k = TRUE} is equivalent to matching with \code{k2k}, then supplying stratum membership as an exact matching variable (i.e., in \code{exact}) to another call to \code{matchit()} with \code{method = "nearest"}, \code{distance = "mahalanobis"} and an argument to \code{discard} denoting unmatched units. It is also equivalent to performing nearest neighbor matching supplying coarsened versions of the variables to \code{exact}, except that \code{method = "cem"} automatically coarsens the continuous variables. The \code{estimand} argument supplied with \code{method = "cem"} functions the same way it would in these alternate matching calls, i.e., by determining the "focal" group that controls the order of the matching.
}
\references{
In a manuscript, you don't need to cite another package when using \code{method = "cem"} because the matching is performed completely within \pkg{MatchIt}. For example, a sentence might read:

\emph{Coarsened exact matching was performed using the MatchIt package (Ho, Imai, King, & Stuart, 2011) in R.}

It would be a good idea to cite the following article, which develops the theory behind coarsened exact matching:

Iacus, S. M., King, G., & Porro, G. (2012). Causal Inference without Balance Checking: Coarsened Exact Matching. Political Analysis, 20(1), 1–24. \doi{10.1093/pan/mpr013}

}

\seealso{
\code{\link{matchit}} for a detailed explanation of the inputs and outputs of a call to \code{matchit}.

The \pkg{cem} package, upon which this method is based and which provided the workhorse in previous versions of \pkg{MatchIt}.

\code{\link{method_exact}} for exact matching, which performs exact matching on the covariates without coarsening.
}
\examples{
data("lalonde")

# Coarsened exact matching on age, race, married, and educ with educ
# coarsened into 5 bins and race coarsened into 2 categories,
# grouping "white" and "hispan" together
m.out1 <- matchit(treat ~ age + race + married + educ, data = lalonde,
                  method = "cem", cutpoints = list(educ = 5),
                  grouping = list(race = list(c("white", "hispan"),
                                              "black")))
m.out1
summary(m.out1)

# The same but requesting 1:1 Mahalanobis distance matching with
# the k2k and k2k.method argument. Note the remaining number of units
# is smaller than when retaining the full matched sample.
m.out2 <- matchit(treat ~ age + race + married + educ, data = lalonde,
                  method = "cem", cutpoints = list(educ = 5),
                  grouping = list(race = list(c("white", "hispan"),
                                              "black")),
                  k2k = TRUE, k2k.method = "mahalanobis")
m.out2
summary(m.out2)
}