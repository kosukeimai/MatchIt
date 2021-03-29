\name{match.data}
\alias{match.data}
\alias{get_matches}

\title{
Construct a matched dataset from a \code{matchit} object
}
\description{
\code{match.data} and \code{get_matches} create a data frame with additional variables for the distance measure, matching weights, and subclasses after matching. This dataset can be used to estimate treatment effects after matching or subclassification. \code{get_matches} is most useful after matching with replacement; otherwise, \code{match.data} is more flexible. See Details below for the difference between them.
}
\usage{
match.data(object,
           group = "all",
           distance = "distance",
           weights = "weights",
           subclass = "subclass",
           data = NULL,
           include.s.weights = TRUE,
           drop.unmatched = TRUE)

get_matches(object,
            distance = "distance",
            weights = "weights",
            subclass = "subclass",
            id = "id",
            data = NULL,
            include.s.weights = TRUE)
}

\arguments{
  \item{object}{
a \code{matchit} object; the output of a call to \fun{matchit}.
}
  \item{group}{
which group should comprise the matched dataset: \code{"all"} for all units, \code{"treated"} for just treated units, or \code{"control"} for just control units. Default is \code{"all"}.
}
  \item{distance}{
a string containing the name that should be given to the variable containing the distance measure in the data frame output. Default is \code{"distance"}, but \code{"prop.score"} or similar might be a good alternative if propensity scores were used in matching. Ignored if a distance measure was not supplied or estimated in the call to \code{matchit}.
}
  \item{weights}{
a string containing the name that should be given to the variable containing the matching weights in the data frame output. Default is \code{"weights"}.
}
  \item{subclass}{
a string containing the name that should be given to the variable containing the subclasses or matched pair membership in the data frame output. Default is \code{"subclass"}.
}
  \item{id}{
a string containing the name that should be given to the variable containing the unit IDs in the data frame output. Default is \code{"id"}. Only used with \code{get_matches}; for \code{match.data}, the units IDs are stored in the row names of the returned data frame.
}
  \item{data}{
a data frame containing the original dataset to which the computed output variables (\code{distance}, \code{weights}, and/or \code{sublcass}) should be appended. If empty, \code{match.data} will attempt to find the dataset using the environment of the \code{matchit} object, which can make this unreliable if \code{match.data} is used in a fresh R session or environment different from the original calling environment (e.g., inside a function) or if the original dataset changed between calling \code{matchit} and \code{match.data}. It is always safest to supply a data frame, which should have as many rows as and be in the same order as the data in the original call to \code{matchit()}. The same goes for \code{get_matches}, which calls \code{match.data} internally.
}
  \item{include.s.weights}{
\code{logical}; whether to multiply the estimated weights by the sampling weights supplied to \code{matchit()}, if any. Default is \code{TRUE}. If \code{FALSE}, the weights in the \code{match.data} or \code{get_matches} output should be multiplied by the sampling weights before being supplied to the function estimating the treatment effect in the matched data.
}
  \item{drop.unmatched}{
\code{logical}; whether the returned data frame should contain all units (\code{FALSE}) or only units that were matched (i.e., have a matching weight greater than zero) (\code{TRUE}). Default is \code{TRUE} to drop unmatched units.
}
}
\details{
\code{match.data} creates a dataset with one row per unit. It will be identical to the dataset supplied except that several new columns will be added containing information related to the matching. When \code{drop.unmatched = TRUE}, the default, units with weights of zero, which are those units that were discarded by common support or the caliper or were simply not matched, will be dropped from the dataset, leaving only the subset of matched units. The idea is for the output of \code{match.data} to be used as the dataset input in calls to \code{glm} or similar to estimate treatment effects in the matched sample. It is important to include the weights in the estimation of the effect and its standard error. The subclass column, when created, contains par or subclass membership and should be used to estimate the effect and its standard error. Subclasses will only be included if there is a \code{subclass} component in the \code{matchit} object, which does not occur with matching with replacement, in which case \code{get_matches} should be used. See \code{vignette("estimating-effects")} for information on how to use \code{match.data} output to estimate effects.

\code{get_matches} is similar to \code{match.data}; the primary difference occurs when matching is performed with replacement, i.e., when units do not belong to a single matched pair. In this case, the output of \code{get_matches} will be a dataset that contains one row per unit for each pair they are a part of. For example, if matching was performed with replacement and a control unit was matched to two treated units, that control unit will have two rows in the output dataset, one for each pair it is a part of. Weights are computed for each row, and are equal to the inverse of the number of control units in each control unit's subclass. Unmatched units are dropped. An additional column with unit IDs will be created (named using the \code{id} argument) to identify when the same unit is present in multiple rows. This dataset structure allows for the inclusion of both subclass membership and repeated use of units, unlike the output of \code{match.data}, which lacks subclass membership when matching is done with replacement. A \code{match.matrix} component of the \code{matchit} object must be present to use \code{get_matches}; in some forms of matching, it is absent, in which case \code{match.data} should be used instead. See \code{vignette("estimating-effects")} for information on how to use \code{get_matches} output to estimate effects after matching with replacement.
}
\value{
A data frame containing the data supplied in the \code{data} argument or in the original call to \code{matchit} with the computed output variables appended as additional columns, named according the arguments above. For \code{match.data}, the \code{group} and \code{drop.unmatched} arguments control whether only subsets of the data are returned. See Details above for how \code{match.data} and \code{get_matches} differ. Note that \code{get_matches} sorts the data by subclass and treatment status, unlike \code{match.data}, which uses the order of the data.

The returned data frame will contain the variables in the original data set or dataset supplied to \code{data}, and the following columns:

\item{distance}{The propensity score, if estimated or supplied to the \code{distance} argument in \code{matchit()}.}
\item{weights}{The computed matching weights. These must be used in effect estimation to correctly incorporate the matching.}
\item{subclass}{Matching strata membership. Units with the same value are in the same stratum.}
\item{id}{The ID of each unit, corresponding to the row names in the original data or dataset supplied to \code{data}. Only included in \code{get_matches} output. This column can be used to identify which rows belong to the same unit since the same unit may appear multiple times if reused in matching with replacement.}

These columns will take on the name supplied to the corresponding arguments in the call to \code{match.data} or \code{get_matches}. See Examples for an example of rename the \code{distance} column to \code{"prop.score"}.

If \code{data} or the original dataset supplied to \code{matchit} was a \code{data.table} or \code{tbl}, the \code{match.data} output will have the same class, but the \code{get_matches} output will always be a base R \code{data.frame}.

In addition to their base class (e.g., \code{data.frame} or \code{tbl}), returned objects have the class \code{matchdata} or \code{getmatches}. This class is important when using \code{\link[=rbind.matchdata]{rbind()}} to append matched datasets.
}

\seealso{
\fun{matchit}

\fun{rbind.matchdata}

\code{vignette("estimating-effects")} for uses of \code{match.data()} and \code{get_matches()} in estimating treatment effects.
}
\examples{
data("lalonde")

# 4:1 matching w/replacement
m.out1 <- matchit(treat ~ age + educ + married +
                    race + nodegree + re74 + re75,
                  data = lalonde, replace = TRUE,
                  caliper = .05, ratio = 4)

m.data1 <- match.data(m.out1, data = lalonde,
                      distance = "prop.score")
dim(m.data1) #one row per matched unit
head(m.data1, 10)

g.matches1 <- get_matches(m.out1, data = lalonde,
                          distance = "prop.score")
dim(g.matches1) #multiple rows per matched unit
head(g.matches1, 10)
}