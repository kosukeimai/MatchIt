\name{method_cem}
\alias{method_cem}

\title{Coarsened Exact Matching}

\description{
In \fun{matchit}, setting \code{method = "cem"} performs coarsened exact matching. With coarsened exact matching, covariates are coarsened into bins, and a complete cross of the coarsened covariates is used to form subclasses defined by each combination of the coarsened covariate levels. Any subclass that doesn't contain both treated and control units is discarded, leaving only subclasses containing treatment and control units that are exactly equal on the coarsened covariates. The coarsening process can be controlled by an algorithm or by manually specifying cutpoints and groupings. The benefits of coarsened exact matching are that the tradeoff between exact matching and approximate balancing can be managed to prevent discarding too many units, which can otherwise occur with exact matching.

This page details the allowable arguments with \code{method = "cem"}. See \fun{matchit} for an explanation of what each argument means in a general context and how it can be specified.

Below is how \code{matchit} is used for coarsened exact matching:
\preformatted{
matchit(formula,
        data = NULL,
        method = "cem",
        estimand = "ATT",
        s.weights = NULL,
        verbose = FALSE,
        ...)
}
}
\arguments{
   \item{formula}{
a two-sided \fun{formula} object containing the treatment and covariates to be used in creating the subclasses defined by a full cross of the coarsened covariate levels.
}
  \item{data}{
a data frame containing the variables named in \code{formula}. If not found in \code{data}, the variables will be sought in the environment.
}
  \item{method}{
set here to \code{"cem"}.
}
  \item{estimand}{
a string containing the desired estimand. Allowable options include \code{"ATT"}, \code{"ATC"}, and \code{"ATE"}. The estimand controls how the weights are computed; see the Computing Weights section at \fun{matchit} for details. When \code{k2k = TRUE} (see below), \code{estimand} also controls how the matching is done.
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
  \item{\code{grouping}}{
     a named list with an (optional) entry for each categorical variable to be matched on. Each element should itself be a list, and each entry of the sublist should be a vector containing levels of the variable that should be combined to form a single level. Any categorical variables not included in \code{grouping} will remain as they are in the data, which means exact matching, with no coarsening, will take place on these variables. See Details.
     }
     \item{\code{cutpoints}}{
        a named list with an (optional) entry for each numeric variable to be matched on. Each element describes a way of coarsening the corresponding variable. They can be a vector of cutpoints that demarcate bins, a single number giving the number of bins, or a string corresponding to a method of computing the number of bins. Allowable strings include \code{"sturges"}, \code{"scott"}, and \code{"fd"}, which use the functions \pkgfun2{grDevices}{nclass}{nclass.Sturges}, \pkgfun2{grDevices}{nclass}{nclass.scott}, and \pkgfun2{grDevices}{nclass}{nclass.FD}, respectively. The default is \code{"sturges"} for variables that are not listed or if no argument is supplied. Can also be a single value to be applied to all numeric variables. See Details.
     }
     \item{\code{k2k}}{
     code{logical}; whether 1:1 matching should occur within the matched strata. If \code{TRUE} nearest neighbor matching without replacement will take place within each stratum, and any unmatched units will be dropped (e.g., if there are more treated than control units in the stratum, the treated units without a match will be dropped). The \code{k2k.method} argument controls how the distance between units is calculated.
     }
     \item{\code{k2k.method}}{
     \code{character}; how the distance between units should be calculated if \code{k2k = TRUE}. Allowable arguments include \code{NULL} (for random matching), \code{"mahalanobis"} (for Mahalanobis distance matching), or any allowable argument to \code{method} in \fun{dist}. Matching will take place on scaled versions of the original (non-coarsened) variables. The default is \code{"mahalanobis"}.
     }
     \item{\code{mpower}}{
     if \code{k2k.method = "minkowski"}, the power used in creating the distance. This is passed to the \code{p} argument of \fun{dist}.
     }
  }
}
The arguments \code{distance} (and related arguments), \code{exact}, \code{mahvars}, \code{discard} (and related arguments), \code{replace}, \code{m.order}, \code{caliper} (and related arguments), and \code{ratio} are ignored with a warning.
}
\section{Outputs}{
All outputs described in \fun{matchit} are returned with \code{method = "cem"} except for \code{match.matrix}. When \code{k2k = TRUE}, a \code{match.matrix} component with the matched pairs is also included.
}
\details{
If the coarsening is such that there are no exact matches with the coarsened variables, the \code{grouping} and \code{cutpoints} arguments can be used to modify the matching specification. Reducing the number of cutpoints or grouping some variable values together can make it easier to find matches. See Examples below. Removing variables can also help (but they will likely not be balanced unless highly correlated with the included variables). To take advantage of coarsened exact matching without failing to find any matches, the covariates can be manually coarsened outside of \code{matchit()} and then supplied to the \code{exact} argument in a call to \code{matchit()} with another matching method.

Setting \code{k2k = TRUE} is equivalent to matching with \code{k2k = FALSE} and then supplying stratum membership as an exact matching variable (i.e., in \code{exact}) to another call to \code{matchit()} with \code{method = "nearest"}, \code{distance = "mahalanobis"} and an argument to \code{discard} denoting unmatched units. It is also equivalent to performing nearest neighbor matching supplying coarsened versions of the variables to \code{exact}, except that \code{method = "cem"} automatically coarsens the continuous variables. The \code{estimand} argument supplied with \code{method = "cem"} functions the same way it would in these alternate matching calls, i.e., by determining the "focal" group that controls the order of the matching.

\subsection{Grouping and Cutpoints}{
The \code{grouping} and \code{cutpoints} arguments allow one to fine-tune the coarsening of the covariates. \code{grouping} is used for combining categories of categorical covariates and \code{cutpoints} is used for binning numeric covariates. The values supplied to these arguments should be iteratively changed until a matching solution that balances covariate balance and remaining sample size is obtained. The arguments are described below.

The argument to \code{grouping} must be a list, where each component has the name of a categorical variable, the levels of which are to be combined. Each component must itself be a list; this list contains one or more vectors of levels, where each vector corresponds to the levels that should be combined into a single category. For example, if a variable \code{amount} had levels \code{"none"}, \code{"some"}, and \code{"a lot"}, one could enter \code{grouping = list(amount = list(c("none"), c("some", "a lot")))}, which would group \code{"some"} and \code{"a lot"} into a single category and leave \code{"none"} in its own category. Any levels left out of the list for each variable will be left alone (so \code{c("none")} could have been omitted from the previous code). Note that if a categorical variable does not appear in \code{grouping}, it will not be coarsened, so exact matching will take place on it. \code{grouping} should not be used for numeric variables; use \code{cutpoints}, described below, instead.

The argument to \code{cutpoints} must also be a list, where each component has the name of a numeric variables that is to be binned. (As a shortcut, it can also be a single value that will be applied to all numeric variables). Each component can take one of three forms: a vector of cutpoints that separate the bins, a single number giving the number of bins, or a string corresponding to an algorithm used to compute the number of bins. Any values at a boundary will be placed into the higher bin; e.g., if the cutpoints were \code{(c(0, 5, 10))}, values of 5 would be placed into the same bin as values of 6, 7, 8, or 9, and values of 10 would be placed into a different bin. Internally, values of \code{-Inf} and \code{Inf} are appended to the beginning and end of the range. When given as a single number defining the number of bins, the bin boundaries are the maximum and minimum values of the variable with bin boundaries evenly spaced between them, i.e., not quantiles. A value of 0 will not perform any binning (equivalent to exact matching on the variable), and a value of 1 will remove the variable from the exact matching variables but it will be still used for pair matching when \code{k2k = TRUE}. The allowable strings include \code{"sturges"}, \code{"scott"}, and \code{"fd"}, which use the corresponding binning method, and \code{"q#"} where \code{#} is a number, which splits the variable into \code{#} equally-sized bins (i.e., quantiles).

An example of a way to supply an argument to \code{cutpoints} would be the following:
\preformatted{
cutpoints = list(X1 = 4,
                 X2 = c(1.7, 5.5, 10.2),
                 X3 = "scott",
                 X4 = "q5")
}
This would split \code{X1} into 4 bins, \code{X2} into bins based on the provided boundaries, \code{X3} into a number of bins determined by \pkgfun2{grDevices}{nclass}{nclass.scott}, and \code{X4} into quintiles. All other numeric variables would be split into a number of bins determined by \pkgfun2{grDevices}{nclass}{nclass.Sturges}, the default.
}
}
\note{
This method does not rely on the \emph{cem} package, instead using code written for \emph{MatchIt}, but its design is based on the original \emph{cem} functions. Versions of \emph{MatchIt} prior to 4.1.0 did rely on \emph{cem}, so results may differ between versions. There are a few differences between the ways \emph{MatchIt} and \emph{cem} (and older versions of \emph{MatchIt}) differ in executing coarsened exact matching, described below.
\itemize{
\item{In \emph{MatchIt}, when a single number is supplied to \code{cutpoints}, it describes the number of bins; in \emph{cem}, it describes the number of cutpoints separating bins. The \emph{MatchIt} method is closer to how \fun{hist} processes breaks points to create bins.}
\item{When \code{cutpoints} are used, \code{"ss"} (for Shimazaki-Shinomoto's rule) can be used in \emph{cem} but not in \emph{MatchIt}.}
\item{When \code{k2k = TRUE}, \emph{MatchIt} matches on the original variables (scaled), whereas \emph{cem} matches on the coarsened variables. Because the variables are already exactly matched on the coarsened variables, matching in \emph{cem} is equivalent to random matching within strata.}
\item{When \code{k2k = TRUE}, in \emph{MatchIt} matched units are identified by pair membership, and the original stratum membership prior to 1:1 matching is discarded. In \emph{cem}, pairs are not identified beyond the stratum the members are part of.}
\item{When \code{k2k = TRUE}, \code{k2k.method = "mahalanobis"} can be requested in \emph{MatchIt} but not in \emph{cem}.}
}
}
\references{
In a manuscript, you don't need to cite another package when using \code{method = "cem"} because the matching is performed completely within \emph{MatchIt}. For example, a sentence might read:

\emph{Coarsened exact matching was performed using the MatchIt package (Ho, Imai, King, & Stuart, 2011) in R.}

It would be a good idea to cite the following article, which develops the theory behind coarsened exact matching:

Iacus, S. M., King, G., & Porro, G. (2012). Causal Inference without Balance Checking: Coarsened Exact Matching. Political Analysis, 20(1), 1–24. \doi{10.1093/pan/mpr013}

}

\seealso{
\fun{matchit} for a detailed explanation of the inputs and outputs of a call to \code{matchit}.

The \emph{cem} package, upon which this method is based and which provided the workhorse in previous versions of \emph{MatchIt}.

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
                                              c("black"))))
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