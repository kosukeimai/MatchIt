\name{summary.matchit}
\alias{summary.matchit}
\alias{summary.matchit.subclass}
\alias{print.summary.matchit}
\alias{print.summary.matchit.subclass}
\title{
View a balance summary of a \code{matchit} object
}
\description{
Computes and prints balance statistics for \code{matchit} and \code{matchit.subclass} objects. Balance should be assessed to ensure the matching or subclassification was effective at eliminating treatment group imbalance and should be reported in the write-up of the results of the analysis.
}
\usage{
\method{summary}{matchit}(object, interactions = FALSE,
        addlvariables = NULL, standardize = TRUE,
        data = NULL, pair.dist = TRUE, ...)
\method{summary}{matchit.subclass}(object, interactions = FALSE,
        addlvariables = NULL, standardize = TRUE,
        data = NULL, pair.dist = FALSE, subclass = FALSE, ...)

\method{print}{summary.matchit}(x, digits = max(3, getOption("digits") - 3),
      ...)
\method{print}{summary.matchit.subclass}(x, digits = max(3, getOption("digits") - 3),
      ...)
}
\arguments{
  \item{object}{
a \code{matchit} object; the output of a call to \code{\link{matchit}}.
}
  \item{interactions}{
\code{logical}; whether to compute balance statistics for two-way interactions and squares of covariates. Default is \code{FALSE}.
}
  \item{addlvariables}{
additional variable for which balance statistics are to be computed along with the covariates in the \code{matchit} object. Can be entered in one of three ways: as a data frame of covariates with as many rows as there were units in the original \code{matchit} call, as a string containing the names of variables in \code{data}, or as a right-sided \code{formula} with the additional variables (and possibly their transformations) found in \code{data}, the environment, or the \code{matchit} object. Balance on squares and interactions of the additional variables will be included if \code{interactions = TRUE}.
}
  \item{standardize}{
\code{logical}; whether to compute standardized (\code{TRUE}) or unstandardized (\code{FALSE}) statistics. The standardized statistics are the standardized mean difference and the median, mean, and maximum of the difference in the (weighted) empirical cumulative distribution functions (ECDFs). The unstandardized statistics are the raw mean difference and the median, mean, and maximum of the quantile-quantile (QQ) difference. See Details below. Default is \code{TRUE}.
}
  \item{data}{
a optional data frame containing variables named in \code{addlvariables} if specified as a string or formula.
}
  \item{pair.dist}{
\code{logical}; whether to compute average absolute pair distances. For matching methods that don't include a \code{match.matrix} component in the output (i.e., exact matching, coarsened exact matching, full matching, and subclassification), computing pair differences can take a long time, especially for large datasets and with many covariates. For other methods (i.e., nearest neighbor, optimal, and genetic matching), computation is fairly quick. Default is \code{FALSE} for subclassification and \code{TRUE} otherwise.
}
  \item{subclass}{
after subclassification, whether to display balance for individual subclasses, and, if so, for which ones. Can be \code{TRUE} (display balance for all subclasses), \code{FALSE} (display balance only in aggregate), or the indices (e.g., \code{1:6}) of the specific subclasses for which to display balance. When anything other than \code{FALSE}, aggregate balance statistics will not be displayed. Default is \code{FALSE}.
}
  \item{digits}{
the number of digits to round balance statistics to.
}
  \item{x}{
a \code{summay.matchit} or \code{summary.matchit.subclass} object; the output of a call to \code{summary}.
}
  \item{\dots}{
ignored.
}
}
\details{
\code{summary} computes a balance summary of a \code{matchit} object. This include balance before and after matching or subclassification, as well as the percent improvement in balance. The variables for which balance statistics are computed are those included in the \code{formula}, \code{exact}, and \code{mahvars} arguments to \code{\link{matchit}}, as well as the distance measure if \code{distance} is not \code{"mahalanobis"}. The \code{X} component of the \code{matchit} object is used to supply the covariates.

The standardized mean differences are computed both before and after matching or subclassification as the difference in treatment group means divided by a standardization factor computed in the unmatched (original) sample. The standardization factor depends on the argument supplied to \code{estimand} in \code{matchit}: for \code{"ATT"}, it is the standard deviation in the treated group; for \code{"ATC"}, it is the standard deviation in the control group; for \code{"ATE"}, it is the square root of the average of the variances within each treatment group. The post-matching mean difference is computed with weighted means in the treatment groups using the matching or subclassification weights.

The variance ratio is computed as the ratio of the treatment group variances. Variance ratios are not computed for binary variables because their variance is a function solely of their mean. After matching, weighted variances are computed using the formula used in \code{\link{cov.wt}}. The percent reduction in bias is computed using the log of the variance ratios.

The eCDF difference statistics are computed by creating a (weighted) eCDF for each group and taking the difference between them for each covariate value. The eCDF is a function that outputs the (weighted)  proportion of units with covariate values at or lower than the input value. The maximum eCDF difference is the same thing as the Kolmogorov-Smirnov statistic. The values are bounded at zero and one, with values closer to zero indicating good overlap between the covariate distributions in the treated and control groups. For binary variables, all eCDF differences are equal to the (weighted) difference in proportion and are computed that way.

The QQ difference statistics are computed by creating two samples of the same size by interpolating the values of the larger one. The values are arranged in order for each sample. The QQ difference for each quantile is the difference between the observed covariate values at that quantile between the two groups. The difference is on the scale of the original covariate. Values close to zero indicate good overlap between the covariate distributions in the treated and control groups. A weighted interpolation is used for post-matching QQ differences. For binary variables, all QQ differences are equal to the (weighted) difference in proportion and are computed that way.

The pair distance is the average of the absolute differences of a variable between pairs. For example, if a treated unit was paired with four control units, that set of units would contribute four absolute differences to the average. Within a subclass, each combination of treated and control unit forms a pair that contributes once to the average. The pair distance is described in Stuart and Green (2008) and is the value that is minimized when using optimal (full) matching. When \code{standardize = TRUE}, the standardized versions of the variables are used, where the standardization factor is as described above for the standardized mean differences. Pair distances are not computed in the unmatched sample (because there are no pairs). Because pair distance can take a while to compute, especially with large datasets or for many covariates, setting \code{pair.dist = FALSE} is one way to speed up \code{summary}.

The effective sample size (ESS) is a measure of the size of a hypothetical unweighted sample with roughly the same precision as a weighted sample. When non-uniform matching weights are computed (e.g., as a result of full matching, matching with replacement, or subclassification), the ESS can be used to quantify the potential precision remaining in the matched sample. The ESS will always be less than or equal to the matched sample size, reflecting the loss in precision due to using the weights. With non-uniform weights, it is printed in the sample size table; otherwise, it is removed because it does not contain additional information above the matched sample size.

After subclassification, the aggregate balance statistics are computed using the subclassification weights rather than averaging across subclasses.

All balance statistics (except pair differences) are computed incorporating the sampling weights upplied to \code{matchit}, if any. The unadjusted balance statistics include the sampling weights and the adjusted balance statistics use the matching weights multiplied by the sampling weights.

When printing, \code{NA} values are replaced with periods (\code{.}), and the pair distance column in the unmatched and percent balance improvement components of the output are omitted.
}
\value{
For \code{matchit} objects, a \code{summary.matchit} object, which is a list with the following components:
\itemize{
  \item{\code{call} - the original call to \code{\link{matchit}}}
  \item{\code{nn} - a matrix of the sample sizes in the original (unmatched) and matched samples}
  \item{\code{sum.all} - a matrix of balance statistics for each covariate in the original (unmatched) sample}
  \item{\code{sum.matched} - a matrix of balance statistics for each covariate in the matched sample}
  \item{\code{reduction} - a matrix of the percent reduction in imbalance for each covariate in the matched sample}
}
For \code{match.subclass} objects, a \code{summary.matchit.subclass} object, which is a list as above with the following additional components:
\itemize{
  \item{\code{call} - the original call to \code{\link{matchit}}}
  \item{\code{sum.all} - a matrix of balance statistics for each covariate in the original sample}
  \item{\code{sum.subclass} - if \code{subclass} is not \code{FALSE}, a list of matrices of balance statistics for each subclass}
  \item{\code{sum.across} - a matrix of balance statistics for each covariate computed using the subclassification weights}
  \item{\code{reduction} - a matrix of the percent reduction in imbalance for each covariate in the matched sample}
  \item{\code{qn} - a matrix of sample sizes within each subclass}
  \item{\code{nn} - a matrix of the sample sizes in the original (unmatched) and matched samples}
}
}
\seealso{
\code{\link{summary}} for the generic method; \code{\link{plot.summary.matchit}} for making a Love plot from \code{summary} output.

\pkgfun2{cobalt}{bal.tab.matchit}{bal.tab} in the \pkg{cobalt} package, which also displays balance for \code{matchit} objects.
}
\examples{
data("lalonde")
m.out <- matchit(treat ~ age + educ + married +
                   race + re74, data = lalonde,
                 method = "nearest", exact = ~ married,
                 replace = TRUE)
summary(m.out, interactions = TRUE)

s.out <- matchit(treat ~ age + educ + married +
                   race + nodegree + re74 + re75,
                 data = lalonde, method = "subclass")
summary(s.out, addlvariables = ~log(age) + I(re74==0))
summary(s.out, subclass = TRUE)
}