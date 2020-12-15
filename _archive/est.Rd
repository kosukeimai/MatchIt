\name{est}
\alias{est_effect}
\alias{print.eff_matchit}
\alias{coef.eff_matchit}
\alias{vcov.eff_matchit}
\title{
Estimate treatment effects after matching
}
\description{
After matching, checking balance, and choosing a single matching specification to use to estimate the effect, it comes time to estimate the effect of the treatment on an outcome. \code{est_effect} makes this process straightforward, requiring little additional coding. \code{est_effect} estimates the marginal effect of treatment and its confidence interval.

Despite the ease of estimating effects \code{est_effect} offers, it is important to remember that it should only be used once (or perhaps twice for a sensitivity analysis) per treatment effect estimated. Running \code{est_effect} many times on the same treatment and outcome to estimate effects with different matching specifications or outcome models is unethical and can invalidate results. Although matching in general reduces the dependence of the effect on the specified model, model dependence can still remain.
}
\usage{
est_effect(formula, m, type, measure, nboot, boot.type, robust.type = NULL, conf = .95, ...)

\method{print}{eff_matchit}(x, digits = max(3, getOption("digits") - 3), ...)
\method{coef}{eff_matchit}(object, link = TRUE, ...)
\method{vcov}{eff_matchit}(object, ...)
}
\arguments{
  \item{formula}{
a \code{\link{formula}} with the outcome variable on the left hand side and the treatment variable (and possibly covariates and treatment-covariate interactions) on the right-hand side. Including covariates in the outcome model is generally recommended when possible, as doing so usually improves precision and can reduce bias.
}
  \item{m}{
a \code{matchit}; the result of a call to \code{\link{matchit}}. Running \code{\link{match.data}} on this object (possibly with other arguments passed through \dots) must produce a dataset containing all the variables named in \code{formula}.
}
  \item{type}{
the method used to estimate the standard error and confidence interval. Allowable options include \code{"clusterrobust"} for cluster-robust standard errors (with pair membership as the cluster), \code{"robust"} for robust standard errors that don't account for pairs, \code{"blockboot"} for the block bootstrap, and \code{"boot"} for the full (traditional) bootstrap. Not all methods are compatible with all inputs. See Details. The default after exact matching, coarsened exact matching, and propensity score subclassification is \code{"robust"}; otherwise the default is \code{"clusterrobust"}. Abbreviations are allowed.
}
  \item{measure}{
the effect measure used to quantify the effect of treatment. For continuous outcomes, the only allowed option is \code{"MD"} for the mean difference. For binary outcomes, the allowed options are \code{"RD"} for the risk difference, \code{"RR"} for the risk ratio, and \code{"OR"} for the odds ratio. For survival outcomes, the only allowed option is \code{"HR"} for the hazard ratio. All effect measures are marginal rather than conditional.
}
  \item{nboot}{
when \code{type} is \code{"blockboot"} or \code{"boot"}, the number of bootstrap replications. This value is supplied to the \code{R} argument in \code{\link[boot:boot]{boot}}.
}
  \item{boot.type}{
the method used to construct bootstrap confidence intervals. The allowable options are \code{"norm"}, \code{"basic"}, \code{"perc"}, and \code{"bca"}. This value is supplied to the \code{type} argument in \code{\link[boot:boot.ci]{boot.ci}}. The default is \code{perc}. To use \code{"bca"}, \code{nboot} must be greater than the original sample size when \code{type = "boot"} and must be greater than the number of pairs/strata when \code{type = "blockboot"}.
}
  \item{robust.type}{
when \code{type} is \code{"clusterrobust"} or \code{"robust"}, the type of standard error to estimate. This argument is supplied to the \code{type} argument in \code{\link[sandwich:vcovCL]{vcovCL}} or \code{\link[sandwich:vcovHC]{vcovHC}}, respectively. See their help files for details. Ignored for survival outcomes. When \code{type = "clusterrobust"}, the default is \code{"HC1"} when \code{measure} is \code{"MD"} or \code{"RD"} and \code{"HC0"} when \code{measure} is \code{"RR"} or \code{"OR"}; when \code{type = "robust"}, the default is \code{"HC3"} for all measures.
}
  \item{conf}{
the confidence level used for the confidence intervals. Default is .95. All intervals are two-sided.
}
  \item{...}{
for \code{est_effect}, arguments passed to \code{\link{match.data}} and \code{\link{get_matches}} (e.g., to ensure no variable name conflicts occur or to supply an alternate dataset) and to \code{\link[boot:boot]{boot}} to control running of the bootstrap. Ignored for other functions documented here.
}
  \item{x, object}{
an \code{eff_matchit} object; the output of a call to \code{est_effect}.
}
  \item{digits}{
the number of significant digits to print. This argument is passed to \code{\link{printCoefmat}}.
}
  \item{link}{
whether to return the estimates on the link scale (\code{TRUE}; default) or outcome scale (\code{FALSE}). See Details.
}
}
\details{
\code{est_effect} estimates the marginal effect of treatment on the outcome using the matched dataset extracted from the supplied \code{matchit} object. For continuous and binary outcomes, a version of g-computation is used to estimate the marginal effects and expected potential outcomes under each treatment. For survival outcomes, a Cox proportional hazards model implemented in \code{\link[survival:coxph]{coxph}} is used.

With continuous outcomes, a linear model is used to estimate the treatment effect. \code{formula} is passed to \code{\link{lm}} to fit the outcome model, and then linear contrasts are used to estimate the effect of treatment and the average potential outcomes under each treatment, averaging across the matched sample. When \code{type} is \code{"clusterrobust"} or \code{"robust"}, the corresponding variance-covariance matrix computed by \code{\link[sandwich:vcovCL]{vcovCL}} or \code{\link[sandwich:vcovHC]{vcovHC}} is used to estimate the standard errors, which are used in standard Wald tests and normal-theory confidence intervals. The critical test statistic and the p-value for the test are computed using a T-distribution with the degrees of freedom of the residual of the linear model. When \code{type} is \code{"blockboot"} or \code{"boot"}, the standard errors (whcih are only used in the hypothesis tests) are computed as the standard deviations of the bootstrapped estimates, and the confidence intervals are computed using the method supplied to \code{boot.type}.

With binary outcomes, the model used to estimate the treatment effect depends on the effect measure requested and the presence of covariates in the outcome model. When the outcome model has covariates (i.e., variables other than the treatment, either as main effects or interactions), the effects are estimated using g-computation, with the potential outcomes estimated using logistic regression models. The estimated conditional probabilities under each treatment are averaged across individuals to form average marginal risks, and the effect is computed as the contrast between the average risks based on the desired effect measure \emph{on the link scale}. Bootstrapping must be used to estimate the standard errors and confidence intervals when covariates are included in the outcome models. When the outcome doesn't contain covariates, a generalized linear model is used to estimate the effect on the link scale. Model-based (cluster-)robust standard errors can be computed from the models to be used in constructing confidence intervals and in hypothesis tests, or bootstrapping can be used. For the risk difference, a linear model (i.e., \code{\link{lm}} is used (note the estimate is the same as using a binomial generalized linear model with an identity link but is easier to compute). For the risk ratio, a \code{\link{quasipoisson}} model with a log link is fit using \code{\link{glm}} (note the estimate is the same as using a binomial generalized linear model with a log link but is less likely to have convergence issues). For the odds ratio, a \code{\link{quasibinomial}} model with a logit link is fit using \code{glm}.

With survival outcomes, a Cox proportional hazards model (implemented using \code{\link[survival:coxph]{coxph}) is used to estimate the log hazard ratio. Covariates cannot be included in the outcome model. When \code{type = "robust"}, \code{robust} is set to \code{TRUE} in \code{coxph}, and when \code{type = "clusterrobust"}, pair/stratum membership is additionally supplied to the \code{cluster} argument; the computed variances are used to compute the standard error used in the confidence interval and hypothesis tests. Only the hazard ratio and log hazard ratio (and their corresponding standard errors and confidence intervals) are included in the output of \code{est_effect}.

All effects, standard errors, and confidence intervals are first computed on the link scale of the effect measure. For the mean difference and risk difference, this is just the difference in average potential outcomes. For the risk ratio, odds ratio, and hazard ratio, the log of each measure is computed. To get the effects on the scale of the desired outcome measure, effects and confidence interval limits are exponentiated. All hypothesis tests are performed on the effects on the link scale. Both the effects on the outcome scale and on the link scale are included in the \code{est_effect} output; see the Value section below.

\subsection{Bootstrapping}{
There are two ways to perform bootstrapping with \code{est_effect}: the traditional full bootstrap and the block bootstrap. Both are implemented using the \code{\link[boot:boot]{boot}} function in the \emph{boot} package.

The full bootstrap involves resampling the original dataset and performing matching (and propensity score estimation if included) and effect estimation in each bootstrap sample. The \code{call} component of the \code{matchit} object is used to call \code{matchit} in each sample. This can take a long time with slower matching methods like optimal and full matching. The full bootstrap should probably not be used with genetic matching because of the computational burden. For most matching methods, the full bootstrap produces accurate or conservative confidence intervals (i.e., when the estimator is unbiased, the intervals cover the true effect at or greater than the specified confidence level).

The block bootstrap involves resampling pairs/strata of units that are formed after matching, so matching and propensity score estimation are only performed once. This makes the block bootstrap much faster than the full bootstrap. The block bootstrap standard errors are estimates of cluster-robust standard errors but can be used for any outcome model, making it a good choice in general. Because only pairs are resampled, fewer bootstrap replications are required to use BCa confidence intervals.

The different types of bootstrap confidence intervals (requested using the \code{boot.type} argument) tend to perfom similarly, but the BCa interval (requested with \code{boot.type = "bca"}) is recommended as its estimation error tends to be the lowest of the available methods. To estimate the BCa interval, \code{nboots} must be equal to the size of the original sample when \code{type = "boot"} and equal to the number of pairs/strata when \code{type = "blockboot"}. In general, the greater \code{nboot} is, the more accurate the confidence intervals; more is always better (but slower - this is the tradeoff). Values in the realm of 1000 or greater are recommended. When using \code{boot.type = "perc"}, it is recommended to choose \code{nboot} so that \code{(nboot + 1)*conf} is an integer to prevent interpolation in computing the bounds; usually values of \code{nboot} one less than a multiple of 100 (e.g., 499) work.

To speed up the bootstrapping, parallel computing can be used by setting \code{parallel = "multicore", ncpus = parallel::detectCores()}. See \code{\link[boot:boot]{boot}} for details. (This may not be available on Windows.)
}
\subsection{After matching with replacement}{
When matching with replacement, the way cluster-robust standard errors are computed differs from other methods. Rather than \code{match.data}, \code{\link{get_matches}} is used to generate the dataset used to estimate the effects. This is so that both multipilicity (i.e., the reuse of control units) and pair membership can be accounted for when estimating standard errors with \code{type = "clusterrobust"}. Both subject ID and pair membership are supplied to the \code{cluster} argument of \code{\link[sandwich:vcovCL]{vcovCL}}.
}
}
\value{
for \code{est_effect}, an \code{eff_estimate.matchit} object containing the following three components:
\itemize{
  \item{\code{effects}: a matrix containing the estimated treatment effect requested through \code{measure} and its confidence interval, and, for continuous and binary outcomes, the estimated expected value of the potential outcomes under each treatment level and their confidence intervals.}
  \item{\code{l.effects}: a matrix containing the estimated treatment effect requested through \code{measure} on the "link" scale and its confidence interval, standard error, and test statistics and p-value for the test that the effect differs from zero (on the link scale), and, for continuous and binary outcomes, the estimated expected value of the potential outcomes under each treatment level and their confidence intervals, standard errors, and test results. For the risk ratio, odds ratio, and hazard ratio, the estimates correspond to the log of the effects, and for the mean difference and risk difference, the estmates are identical to those in \code{effects}.}
  \item{\code{vcov}: the estimate variance-covariance matrix of the estimate effects on the link scale. The square root of the diagonal of this matrix is used to compute the standard errors in \code{l.effects}.}
}

The \code{coef} and \code{vcov} methods extract the effect estimates and covariance matrix of the estimates on the link scale. That is, \code{coef} returns \code{est(...)$l.effects[,"Estimate"]} and \code{vcov} returns \code{est(...)$vcov}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Noah Greifer
}
\note{
Although many combinations of standard error estimates and matching methods are supported, many of them have not been statistically validated. See \code{vignette("estimating-effects")} for a summary of the research on estimating effects and their standard errors after various forms of matching.
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data("lalonde")

#NN matching
m.out <- matchit(treat ~ age + educ + married + nodegree +
                   race + re74 + re75, data = lalonde,
                   method = "full")

#Check balance
summary(m.out)

#Estimate the effect w/ cluster-robust SE (default)
est(re78 ~ treat + age + educ + married + nodegree +
       race + re74 + re75, m.out)

#See vignette("estimating-effects") for examples with
#other outcome types
}
