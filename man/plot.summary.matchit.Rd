\name{plot.summary.matchit}
\alias{plot.summary.matchit}

\title{
Generate a Love Plot of Standardized Mean Differences
}
\description{
Generates a Love plot, which is a dot plot with variable names on the y-axis and standardized mean differences on the x-axis. Each point represents the standardized mean difference of the corresponding covariate in the matched or unmatched sample. Love plots are a simple way to display covariate balance before and after matching. The plots are generated using \fun{dotchart} and \fun{points}.
}
\usage{
\method{plot}{summary.matchit}(x, abs = TRUE, var.order = "data",
       threshold = c(.1, .05), position = "bottomright", ...)
}

\arguments{
  \item{x}{
a \code{summary.matchit} object; the output of a call to \fun{summary.matchit}. The \code{standardize} argument must be set to \code{TRUE} (which is the default) in the call to \code{summary}.
}
  \item{abs}{
\code{logical}; whether the standardized mean differences should be displayed in absolute value (\code{TRUE}, default) or not \code{FALSE}.
}
  \item{var.order}{
how the variables should be ordered. Allowable options include \code{"data"}, ordering the variables as they appear in the \code{summary} output; \code{"unmatched"}, ordered the variables based on their standardized mean differences before matching; \code{"matched"}, ordered the variables based on their standardized mean differences after matching; and \code{"alphabetical"}, ordering the variables alphabetically. Default is \code{"data"}. Abbreviations allowed.
}
  \item{threshold}{
numeric values at which to place vertical lines indicating a balance threshold. These can make it easier to see for which variables balance has been achieved given a threshold. Multiple values can be supplied to add multiple lines. When \code{abs = FALSE}, the lines will be displayed on both sides of zero. The lines are drawn with \code{abline} with the linetype (\code{lty}) argument corresponding to the order of the entered variables (see options at \fun{par}). The default is \code{c(.1, .05)} for a solid line (\code{lty = 1}) at .1 and a dashed line (\code{lty = 2}) at .05, indicating acceptable and good balance, respectively. Enter a value as \code{NA} to skip that value of \code{lty} (e.g., \code{c(NA, .05)} to have only a dashed vertical line at .05).
}
  \item{position}{
the position of the legend. Should be one of the allowed keyword options supplied to \code{x} in \fun{legend} (e.g., \code{"right"}, \code{"bottomright"}, etc.). Default is \code{"bottomright"}. Set to \code{NULL} for no legend to be included. Note that the legend will cover up points if you are not careful; setting \code{var.order} appropriately can help in avoiding this.
}
  \item{\dots}{
ignored.
}
}
\details{
For matching methods other than subclassification, \code{plot.summary.matchit} uses \code{x$sum.all[,"Std. Mean Diff."]} and \code{x$sum.matched[,"Std. Mean Diff."]} as the x-axis values. For subclassification, in addition to points for the unadjusted and aggregate subclass balance, numerals representing balance in individual subclasses are plotted if \code{subclass = TRUE} in the call to \code{summary}. Aggregate subclass standardized mean differences are taken from \code{x$sum.across[,"Std. Mean Diff."]} and the subclass-specific mean differences are taken from \code{x$sum.subclass}.
}
\value{
A plot is displayed, and \code{x} is invisibly returned.
}
\author{
Noah Greifer
}
\seealso{
\fun{summary.matchit}, \fun{dotchart}

\pkgfun{cobalt}{love.plot} is a more flexible and sophisticated function to make Love plots and is also natively compatible with \code{matchit} objects.
}
\examples{
data("lalonde")
m.out <- matchit(treat ~ age + educ + married +
                   race + re74, data = lalonde,
                 method = "nearest")
plot(summary(m.out, interactions = TRUE),
     var.order = "unmatched")

s.out <- matchit(treat ~ age + educ + married +
                   race + nodegree + re74 + re75,
                 data = lalonde, method = "subclass")
plot(summary(s.out, subclass = TRUE),
     var.order = "unmatched", abs = FALSE)
}
