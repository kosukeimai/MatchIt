\name{method_cardinality}
\alias{method_cardinality}
\title{
Cardinality Matching
}
\description{
In \fun{matchit}, setting \code{method = "cardinality"} performs cardinality matching and other forms of matching using mixed integer programming. Rather than forming pairs, cardinality matching selects the largest subset of units that satisfies user-supplied balance constraints on mean differences. One of several available optimization programs can be used to solve the mixed integer program. The default is the GLPK library as implemented in the \code{Rglpk} package, but performance can be dramatically improved using Gurobi and the \code{gurobi} package, for which there is a free academic license.

This page details the allowable arguments with \code{method = "cardinality"}. See \fun{matchit} for an explanation of what each argument means in a general context and how it can be specified.

Below is how \code{matchit()} is used for nearest neighbor matching:
\preformatted{
matchit(formula,
        data = NULL,
        method = "cardinality",
        estimand = "ATT",
        exact = NULL,
        discard = "none",
        s.weights = NULL,
        ratio = 1,
        verbose = FALSE, ...)
}
}
\arguments{
   \item{formula}{
a two-sided \fun{formula} object containing the treatment and covariates to be balanced.
}
  \item{data}{
a data frame containing the variables named in \code{formula}. If not found in \code{data}, the variables will be sought in the environment.
}
  \item{method}{
set here to \code{"cardinality"}.
}
  \item{estimand}{
a string containing the desired estimand. Allowable options include \code{"ATT"}, \code{"ATC"}, and \code{"ATE"}. See Details.
}
  \item{exact}{
for which variables exact matching should take place. Seperate optimization will occur within each subgroup of the exact matching variables.
}
  \item{discard}{
a string containing a method for discarding units outside a region of common support.
}
  \item{s.weights}{
the variable containing sampling weights to be incorporated into the optimization. The balance constraints refer to the product of the sampling weights and the matching weights, and the sum of the product of the sampling and matching weights will be maximized.
}
  \item{ratio}{
the desired ratio of control to treated units. Can be set to \code{NA} to maximize sample size without concern for this ratio. See \code{Details}.
}
  \item{verbose}{
\code{logical}; whether information about the matching process should be printed to the console. ###EDIT
}
  \item{\dots}{
additional arguments that control the matching specification:
  \describe{
     \item{\code{tols}}{
      \code{numeric}; a vector of imbalance tolerances for mean differences, one for each covariate in \code{formula}. If only one value is supplied, it is applied to all. See \code{std.tols} below. Default is \code{.05} for standardized mean differences of at most .05 for all covariates between the treatment groups in the matched sample.
     }
     \item{\code{std.tols}}{
        \code{logical}; whether each entry in \code{tols} corresponds to a raw or standardized mean difference. If only one value is supplied, it is applied to all. Default is \code{TRUE} for standardized mean differences. The standardization factor is the pooled standard deviation when \code{estimand = "ATE"}, the standard deviation of the treated group when \code{estimand = "ATT"}, and the standard deviation of the control group when \code{estimand = "ATC"} (the same as used in \fun{summary.matchit}).
     }
     \item{\code{solver}}{
        the name of solver to use to solve the otpimziation problem. Available options include \code{"glpk"}, \code{"symphony"}, and \code{"gurobi"} for GLPK (implemented in the \code{Rglpk} package), SYMPHONY (implemented in the \code{Rsymphony} package), and Gurobi (implemented in the \code{"gurobi"} package). The differences between them are in speed and solving abaility. GLPK (the default) is the easiest to install, but Gurobi is recommended as it consistently outperforms other solvers and can find solutions even when others can't, and in less time. Gurobi is proprietary but can be used with a free trial or academic license.
     }
     \item{\code{time}}{
        the maximum amount of time before the optimization routine aborts, in seconds. Default is 120 (2 minutes). For large problems, this should be set much higher.
     }
  }
}
}
\section{Outputs}{
Most outputs described in \fun{matchit} are returned with \code{method = "cardinality"}. The \code{match.matrix} and \code{subclass} components are omitted because no pairing or subclassification is done.
}
\details{
\subsection{Cardinality and Template Matching}{
Two types of matching are available with \code{method = "cardinality"}: cardinality matching and template matching.

Cardinality matching finds the largest matched set that satisfies the balance constraints between treatment groups, with the additional constraint that the ratio of the number of matched control to matched treated units is equal to \code{ratio} (1 by default), mimicking k:1 matching. When not all treated units are included in the matched set, the estimand no longer corresponds to the ATT, so cardinality matching should be avoided if retaining the ATT is desired. To request cardinality matching, \code{estimand} should be set to \code{"ATT"} or \code{"ATC"} and \code{ratio} should be set to a positive integer. 1:1 cardinality matching is the default method when no arguments are specified.

Template matching finds the largest matched set that satisfies balance constraints between each treatment group and a specified target sample. When \code{estimand = "ATT"}, it will find the largest subset of the control units that satisfies the balance constraints with respect to the treated group, which is left intact. When \code{estimand = "ATE"}, it will find the largest subsets of the treated group and of the control group that are balanced to the overall sample. To request template matching for the ATT, \code{estimand} should be set to \code{"ATT"} and \code{"ratio"} to \code{NA}. To request template matching for the ATE, \code{estimand} should be set to \code{"ATE"} and \code{ratio} can be set either to \code{NA} to maximize the size of each sample independently or to a positive integer to ensure that the ratio of matched control units to matched treated treats is fixed, mimicking k:1 matching. Unlike cardinality matching, template matching retains the requested estimand if a solution is found.

Neither method involves creating pairs in the matched set, but it is possible to perform an additional round of pairing within the matched sample after cardinality matching or template matching for the ATE with a fixed sample size ratio. See Examples for an example of optimal pair matching after cardinality matching. The balance will not change, but additional precision and robustness can be gained by forming the pairs.

The weights are scaled so that the sum of the weights in each group is equal to the number of matched units in the smaller group when cardinality matching or template matching for the ATE, and scaled so that the sum of the weights in the control group is equal to the number of treated units when template matching for the ATT. When the sample sizes of the matched groups is the same (i.e., when \code{ratio = 1}), no scaling is done. Robust standard errors should be used in effect estimation after cardinality or template matching (and cluster-robust standard errors if additional pairing is done in the matched sample). See \code{vignette("estimating-effects")} for more information.
}

\subsection{Specifying Balance Constraints}{
The balance constraints are on the (standardized) mean differences between the matched treatment groups for each covariate. Balance constraints should be set by supplying arguments to \code{tols} and \code{std.tols}. For example, setting \code{tols = .1} and \code{std.tols = TRUE} requests that all the mean differences in the matched sample should be within .1 standard deviations for each covariate. Different tolerances can be set for difference variables; it might be beneficial to constrain the mean differences for highly prognostic covariates more tightly than other variables. For example, one could specify \code{tols = c(.001, .05), std.tols = c(TRUE, FALSE)} to request that the standadized mean difference for the first covariate is less than .001 and the raw mean difference for the second covariate is less than .05. The values should be specified in the order they appear in \code{formula}, except when interactions are present. One can run the following code:

\preformatted{MatchIt:::get_assign(model.matrix(~X1*X2 + X3, data = data))[-1]}

this will output a vector of numbers and the variable to which each number corresponds; the first entry in \code{tols} corresponds to the variable labeled 1, the second to the variable labeled 2, etc.
}

\subsection{Dealing with Errors and Warnings}{
When the optimization cannot be solved at all, or at least within the time frame specified in the argument ot \code{time}, an error or warning will appear. Unfortunately, it is hard to know exactly the cause of the failure and what measures should be taken to rectify it.

A warning that says \code{"The optimizer failed to find an optimal solution in the time alotted. The returned solution may not be optimal."} usually means that an optimal solution may be possible to find with more time, in which case \code{time} should be increased or a faster solver should be used. Even with this warning, a potentially useable solution will be returned, so don't automatically take it to mean the optimization failed. Sometimes, when there are multiple solutions with the same resulting sample size, the optimziers will stall at one of them, not thinking it has found the optimum. The result should be checked to see if it can be used as the solution.

An error that says \code{"The optimization problem may be infeasible."} usually means that there is a issue with the optimization problem, i.e., that there is no possible way to satisfy the constraints. To rectify this, one can try relaxing the constraints by increasing the value of \code{tols} or use another solver. Sometimes Gurobi can solve problems that the other solvers cannot.
}

}
\references{
In a manuscript, you should reference the solver used in the optimization. For example, a sentence might read:

\emph{Cardinality matching was performed using the MatchIt package (Ho, Imai, King, & Stuart, 2011) in R with the optimization performed by GLPK.}

See \code{vignette("matching-methods")} for more literature on cardinality matching.
}

\seealso{
\fun{matchit} for a detailed explanation of the inputs and outputs of a call to \code{matchit()}.

\CRANpkg{designmatch}, which performs cardinality and template matching with many more options and more flexibility. The implementations of cardinality matching differ between \code{MatchIt} and \code{designmatch}, so their results might differ.

\CRANpkg{optweight}, which offers similar functionality but in the context of weighting rather than matching.
}
\examples{
data("lalonde")

#Choose your solver; "gurobi" is best, "glpk" is easiest to install
solver = "glpk"

# 1:1 cardinality matching
m.out1 <- matchit(treat ~ age + educ + re74 + re75,
                  data = lalonde, method = "cardinality",
                  estimand = "ATT", ratio = 1,
                  tols = .1, solver = solver)
m.out1
summary(m.out1)

# Template matching for the ATT
m.out2 <- matchit(treat ~ age + educ + re74 + re75,
                  data = lalonde, method = "cardinality",
                  estimand = "ATT", ratio = NA,
                  tols = .1, solver = solver)
m.out2
summary(m.out2, un = FALSE)

# Template matching for the ATE with a constraint on
# sample size
m.out3 <- matchit(treat ~ age + educ + re74 + re75,
                  data = lalonde, method = "cardinality",
                  estimand = "ATE", ratio = 1,
                  tols = .1, solver = solver)
m.out3
summary(m.out3, un = FALSE)

# Pairing after 1:1 cardinality matching:
m.out4 <- matchit(treat ~ age + educ + re74 + re75,
                  data = lalonde, method = "nearest",
                  distance = "mahalanobis",
                  discard = m.out1$weights == 0)

# Note that balance doesn't change but pair distances
# are lower for the paired-upon variables
summary(m.out4, un = FALSE)

# In these examples, a high tol was used and
# few covariate matched on in order to not take too long;
# with real data, tols should be much lower and more
# covariates included if possible.
}
