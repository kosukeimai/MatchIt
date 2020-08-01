\name{plot.matchit}
\alias{plot.matchit}
\alias{plot.matchit.subclass}
\title{
Generate Balance Plots after Matching and Subclassification
}
\description{
Generates plots displaying distributional balance and overlap on covariates and propensity scores before and after matching and subclassification. For displaying balance solely on covariate standardized mean differences, see \code{\link{plot.summary.matchit}}. The plots here can be used to assess to what degree covariate and propensity score distributions are balanced and how weighting and discarding affect the distribution of propensity scores.
}
\usage{
\method{plot}{matchit}(x, type = "qq", interactive = TRUE,
     which.xs = NULL, ...)

\method{plot}{matchit.subclass}(x, type = "qq", interactive = TRUE,
     which.xs = NULL, subclass = NULL, ...)
}
\arguments{
  \item{x}{
a \code{matchit} object; the output of a call to \code{\link{matchit}}.
}
  \item{type}{
the type of plot to display. Options include \code{"qq"}, \code{"jitter"}, and \code{"histogram"}. See Details. Default is \code{"qq"}. Abbreviatiosn allowed.
}
  \item{interactive}{
\code{logical}; whether the graphs should be displayed ina an interactive way. Only applies for \code{type = "qq"} and \code{type = "jitter"}. See Details.
}
  \item{which.xs}{
with \code{type = "qq"}, for which covariate(s) QQ plots should be displayed. Factor variables should be named by the original variable name rather than the names of individual dummy variables created after expansion with \code{model.matrix}.
}
  \item{subclass}{
with subclassification and \code{type = "qq"}, for which subclass QQ plots should be displayed. If \code{NULL}, QQ plots are displayed for the full sample weighted by the subclassification weights. Ignored when \code{interactive = TRUE}.
}
  \item{\dots}{
arguments passed to \code{\link{plot}} to contro lthe appearance of the plot. Not all options are accepted.
}
}
\details{
\code{plot.matchit} makes one of three different plots dependoing on the argument supplied to \code{type}.

With \code{type = "qq"}, quantile-quantile (QQ) plots are created for each covariate before and after matching. The plots involve interpolating points in the smaller group based on the weighted quantiles of the other group. When points are approximately on the 45-degree line, the distributions in the treatment and control groups are approximately equal. Major deviations indicate departures from distributional balance. With variable with fewer than 5 unique values, points are jittered to more easily visualize counts. When \code{interactive = TRUE}, plots for three variables will be displayed at a time, and the prompt in the console allows you to move on to the next set of variables. When \code{interactive = FALSE}, multiple pages are plotted at the same time, but only the last few variables will be visible in the displayed plot. To see only a few specific variables at a time, use the \code{which.xs} argument to display QQ plots for just those variables. If fewer than three (after expanding factors into their dummies), \code{interactive} is ignored.

With \code{type = "jitter"}, a jitter plot is displayed for distance values before and after matching. This method requires a distance variable (e.g., a propensity score) to have been estimated or supplied in the call to \code{matchit}. The plot displays individuals values for matched and unmatched tretament and control units arranged horizontally by their propensity scores. Points are jitter so counts are easier to see. The size of the points increases when they receive higher weights. When \code{interactive = TRUE}, you can click on points in the graph to identify their rownames and indeces to further probe extreme values, for example. With subclassification, vertical lines representing the subclass boundaries are overlay on the plots.

With \code{type = "histogram"}, a histogram of distance values is displayed for the treatment and control groups before and after matching. This method requires a distance variable (e.g., a propensity score) to have been estimated or supplied in the call to \code{matchit}. With subclassification, vertical lines representing the subclass boundaries are overlay on the plots.
}
\value{
A plot is displayed, and \code{x} is invisibly returned.
}
\note{
Sometimes, bugs in the plotting functions can cause strange layout or size issues. Run \code{\link{frame}()} or \code{\link{dev.off}()} can be used to reset the plotting pane (note the latter will delete any plots in the plot history).
}

\seealso{
\code{\link{summary.matchit}} for numerical summaries of balance, including those that rely on the QQ plots.

\code{\link{plot.summary.matchit}} for plotting standardized mean differences in a Love plot.

\code{\link[cobalt]{bal.plot}} in the \pkg{cobalt} package for displaying distributional balance in several other ways that are more easily customizable and rely on \pkg{ggplot2} objects. \pkg{cobalt} functions natively support \code{matchit} objects.
}
\examples{
data("lalonde")
m.out <- matchit(treat ~ age + educ + married +
                   race + re74, data = lalonde,
                 method = "nearest")
plot(m.out, type = "qq", interactive = FALSE,
     which.xs = c("age", "educ", "married"))
plot(m.out, type = "histogram")

s.out <- matchit(treat ~ age + educ + married +
                   race + nodegree + re74 + re75,
                 data = lalonde, method = "subclass")
plot(s.out, type = "jitter", interactive = FALSE)
}
