\name{method_subclass}
\alias{method_subclass}
%%%%%NEEDS EDITING
\title{
Subclassification
}
\description{
In \code{\link{matchit}}, setting \code{method = "subclass"} performs subclassification on the distance measure (i.e., propensity score). Treatment and control units are placed into subclasses based on quantiles of the propensity score in the treated group, in the control gorup, or overall, depending on the desired estimand. Weights are computed based on the proportion of treated units in each subclass. Subclassification implemented here does not rely on any other package.

This page details the allowable arguments with \code{method = "subclass"}. See \code{\link{matchit}} for an explanation of what each argument means in a general context and how it can be specified.
}
\usage{
matchit(formula, data = NULL, method = "genetic",
        distance = "glm", link = "logit",
        distance.options=list(), estimand = "ATT",
        discard = "none", reestimate = FALSE,
        verbose = FALSE, ...)
}
\arguments{
   \item{formula}{
a two-sided \code{\link{formula}} object containing the treatment and covariates to be used in creating the distance measure used in the matching.
}
  \item{data}{
a data frame containing the variables named in \code{formula}. If not found in \code{data}, the variables will be saught in the environment.
}
  \item{method}{
set here to \code{"nearest"}.
}
  \item{distance}{
the distance measure to be used. See \code{\link{distance}} for allowable options. \code{distance = "mahalanobis"} is not allowed.
}
  \item{link}{
when \code{distance} is specified as a string, an additional argument controlling the link function used in estimating the distance measure. See \code{\link{distance}} for allowable options with each option.
}
  \item{distance.options}{
a named list containing additional arguments supplied to the function that estimates the distance measure as determined by the argument to \code{distance}.
}
  \item{estimand}{
the target \code{estimand}. If \code{"ATT"}, the default, subclasses are formed based on quantiles of the distance measure in the treated group; if \code{"ATC"}, subclasses are formed based on quantiles of the distance measure in the control group; if \code{"ATE"}, subclasses are formed based on quantiles of the distance measure in the full sample. \code{estimand} also affects how the subclassification weights are computed; see Details.
}
  \item{discard}{
a string containing a method for discarding units outside a region of common support.
}
  \item{reestimate}{
if \code{discard} is not \code{"none"}, whether to re-estimate the propensity score in the remaining sample prior to subclassification.
}
  \item{verbose}{
\code{logical}; whether information about the matching process should be printed to the console.
}
  \item{\dots}{
additional arguments that control the subclassification:
  \describe{
     \item{\code{subclass}}{
        either the number of subclasses desired or a vector of quantiles used to divide the distance measures into subclasses. Default is 6.
     }
     \item{\code{min.n}}{
     the minimum number of units of each treatment group that are to be assigned each subclass. If the distance measure is divided in such as way that fewer than \code{min.n} units of a treatment group are assigned a given subclass, units from other subclasses will be reassigned to fill the defficient subclass. Default is 1.
     }
  }
}
}
\section{Outputs}{
All outputs described in \code{\link{matchit}} are returned with \code{method = "subclass"} except that \code{match.matrix} is excluded and an additional component, \code{q.cut}, is included, containing a vector of the distance measure cutpoints used to define the subclasses. Note that when \code{min.n > 0}, the subclass assignments may not strictly obey the quantiles listed in \code{q.cut}.
}
\details{
Subclassification weights are formed in the following way. A new subclass "propensity score" (\code{p}) is computed as the proportion of units in each subclass that are in the treated group, and all units in that subclass are assigned that propensity score. Weights are computed based on the argument to \code{estimand} using the standard formulas for inverse probability weights: for the ATT, weights are 1 for the treated units and \code{p/(1-p)} for the control units; for the ATC, weights are \code{(1-p)/p} for the treated units and 1 for the control units; for the ATE, weights are \code{1/p} for the treated units and \code{1/(1-p)} for the control units. In each treatment group, weights are divided by the sum of the weights in that treatment group to make the weights sum to the number of units in that treatment group.

After subclassification, effect estimates can be computed seprately in the subclasses and combined, or a single marginal effect can be estimated by using the weights in the full sample. When using the weights, the method is sometimes referred to as marginal mean weighting through stratification (MMWS; Hong, 2010) or fine stratification weighting (Desai et al., 2017). The weights can be interpreted just like inverse probability weights.

Changing \code{min.n} can change the quality of the weights. Generally, a low \code{min.w} will yield better balance because subclasses only contain units with relatively similar distance values, but may yield higher variance because extreme weights can occur due to there being few of a treatment group in some subclasses.

Note that subclassification weights can also be estimated using \pkg{WeightIt}, which provides some additional methods for estimating propensity scores. Where propensity score-estimation methods overlap, both methods will yield the same weights.

}
\references{
In a manuscript, you don't need to cite another package when using \code{method = "nearest"} because the matching is performed completely within \pkg{MatchIt}. For example, a sentence might read:

\emph{Propensity score subclassification was performed using the MatchIt package (Ho, Imai, King, & Stuart, 2011) in R.}

It may be a good idea to cite Hong (2010) or Desai et al. (2017) if the treatment effect is estimated using the subclassification weights.

Desai, R. J., Rothman, K. J., Bateman, B. . T., Hernandez-Diaz, S., & Huybrechts, K. F. (2017). A Propensity-score-based Fine Stratification Approach for Confounding Adjustment When Exposure Is Infrequent: Epidemiology, 28(2), 249–257. \doi{10.1097/EDE.0000000000000595}

Hong, G. (2010). Marginal mean weighting through stratification: Adjustment for selection bias in multilevel data. Journal of Educational and Behavioral Statistics, 35(5), 499–531. \doi{10.3102/1076998609359785}

}

\seealso{
\code{\link{matchit}} for a detailed explanation of the inputs and outputs of a call to \code{matchit}.

\code{\link{method_full}} for optimal full matching, which is similar to subclassification except that the number of subclasses and subclass membership are chosen to optimize the within-subclass disiatnce.
}
\examples{
data("lalonde")

# PS subclassification for the ATT with 7 subclasses
s.out1 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75, data = lalonde,
                  method = "subclass", subclass = 7)
summary(s.out1, subclass = TRUE)

# PS subclassification for the ATE with 10 subclasses
# and at least 2 units in each group per subclass
s.out2 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75, data = lalonde,
                  method = "subclass", subclass = 10
                  estimand = "ATE", min.n = 2)
summary(s.out2)
}