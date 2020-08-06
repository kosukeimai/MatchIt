\name{method_cem}
\alias{method_cem}

\title{Coarsened Exact Matching}

\description{
In \code{\link{matchit}}, setting \code{method = "cem"} performs coarsened exact matching. With coarsened exact matching, covariates are coaursened into bins, and a complete cross of the coarsened covariates is used to form subclasses defined by each combination of the coarsened covariate levels. Any subclass that doesn't contain both treated and control units is discarded, leaving only subclasses containing treatment and control units that are exactly equal on the coarsened covariates. The coarseneing process can be controlled by an algorithm or by manually specifying cutpoints and groupings. The benefits of coarsened exact matching are that the tradeoff between exact matching and approximate balancing can be managed to prevent discarding too many units, which can otherwise occur with exact matching.

This page details the allowable arguments with \code{method = "cem"}. See \code{\link{matchit}} for an explanation of what each argument means in a general context and how it can be specified.

Below is how \code{matchit} is used for coarsened exact matching:
\preformatted{
matchit(formula, data = NULL, method = "cem",
        estimand = "ATT", verbose = FALSE, ...)
}
}
\arguments{
   \item{formula}{
a two-sided \code{\link{formula}} object containing the treatment and covariates to be used in creating the subclasses defined by a full cross of the coarsened covariate levels.
}
  \item{data}{
a data frame containing the variables named in \code{formula}. If not found in \code{data}, the variables will be saught in the environment.
}
  \item{method}{
set here to \code{"cem"}.
}
  \item{estimand}{
a string containing the desired estimand. Allowable options include \code{"ATT"}, \code{"ATC"}, and \code{"ATE"}. The estimand controls how the weights are computed; see the Computing Weights section at \code{\link{matchit}} for details.
}
  \item{verbose}{
\code{logical}; whether information about the matching process should be printed to the console.
}
  \item{\dots}{
additional arguments passed to \code{\link[cem]{cem}} in \pkg{cem}. Potentially useful options include \code{cutpoints}, \code{grouping}, \code{k2k}, and \code{method}. To supply and argument to \code{method} in \code{cem}, supply an argument to \code{k2k.method} in \code{matchit} in order to avoid a name conflict with \code{method} in \code{matchit}.
}
}
\section{Outputs}{
All outputs described in \code{\link{matchit}} are returned with \code{method = "cem"} except for \code{match.matrix}. This is because matching strata are not indexed by treated units as they are in some other forms of matching.
}
\references{
In a manuscript, be sure to cite the following paper if using \code{matchit} with \code{method = "cem"}:

Iacus, S., King, G., & Porro, G. (2009). cem: Software for Coarsened Exact Matching. Journal of Statistical Software, 30(1), 1–27. \doi{10.18637/jss.v030.i09}

For example, a sentence might read:

\emph{Coarsened matching was performed using the MatchIt package (Ho, Imai, King, & Stuart, 2011) in R, which calls functions from the cem package (Iacus, King, & Porro, 2009).}

Theory is also developed in the following article:

Iacus, S. M., King, G., & Porro, G. (2012). Causal Inference without Balance Checking: Coarsened Exact Matching. Political Analysis, 20(1), 1–24. \doi{10.1093/pan/mpr013}

}

\seealso{
\code{\link{matchit}} for a detailed explanation of the inputs and outputs of a call to \code{matchit}.

\code{\link{method_exact}} for exact matching, which performs exact matching on the covariates without coarsening.
}
\examples{
data("lalonde")

# Coarsened exact matching on age, race, married, and educ with educ
# coarsened into 10 bins
m.out1 <- matchit(treat ~ age + race + married + educ, data = lalonde,
                  method = "cem", cutpoints = list(educ = 9))
m.out1
summary(m.out1)
}