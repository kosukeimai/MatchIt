\name{method_exact}
\alias{method_exact}

\title{Exact Matching}
\description{
In \code{\link{matchit}}, setting \code{method = "exact"} performs exact matching. With exact matching, a complete cross of the covariates is used to form subclasses defined by each combination of the covariate levels. Any subclass that doesn't contain both treated and control units is discarded, leaving only subclasses containing treatment and control units that are exactly equal on the included covariates. The benefits of exact matching are that confounding due to the covariates included is completely eliminated, regardless of the functional form of the treatment or outcome models. The problem is that typically many units will be discarded, sometimes dramatically reducing precision and changing the target population of inference. To use exact matching in combination with another matching method (i.e., to exact match on some covariates and some other form of matching on others), use the \code{exact} argument with that method.

This page details the allowable arguments with \code{method = "exact"}. See \code{\link{matchit}} for an explanation of what each argument means in a general context and how it can be specified.

Below is how \code{matchit} is used for exact matching:
\preformatted{
matchit(formula, data = NULL, method = "exact",
        estimand = "ATT", verbose = FALSE, ...)
}
}
\arguments{
   \item{formula}{
a two-sided \code{\link{formula}} object containing the treatment and covariates to be used in creating the subclasses defined by a full cross of the covariate levels.
}
  \item{data}{
a data frame containing the variables named in \code{formula}. If not found in \code{data}, the variables will be sought in the environment.
}
  \item{method}{
set here to \code{"exact"}.
}
  \item{estimand}{
a string containing the desired estimand. Allowable options include \code{"ATT"}, \code{"ATC"}, and \code{"ATE"}. The estimand controls how the weights are computed; see the Computing Weights section at \code{\link{matchit}} for details.
}
  \item{verbose}{
\code{logical}; whether information about the matching process should be printed to the console.
}
  \item{\dots}{
ignored.
}
}
\section{Outputs}{
All outputs described in \code{\link{matchit}} are returned with \code{method = "exact"} except for \code{match.matrix}. This is because matching strata are not indexed by treated units as they are in some other forms of matching.
}
\references{
In a manuscript, you don't need to cite another package when using \code{method = "exact"} because the matching is performed completely within \pkg{MatchIt}. For example, a sentence might read:

\emph{Exact matching was performed using the MatchIt package (Ho, Imai, King, & Stuart, 2011) in R.}
}

\seealso{
\code{\link{matchit}} for a detailed explanation of the inputs and outputs of a call to \code{matchit}. The \code{exact} argument can be used with other methods to perform exact matching in combination with other matching methods.

\code{\link{method_cem}} for coarsened exact matching, which performs exact matching on coarsened versions of the covariates.
}
\examples{
data("lalonde")

# Exact matching on age, race, married, and educ
m.out1 <- matchit(treat ~ age + race + married + educ, data = lalonde,
                  method = "exact")
m.out1
summary(m.out1)
}