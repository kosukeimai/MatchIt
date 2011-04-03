\name{match.data}

\alias{match.data}

\title{Output Matched Data Sets}

\description{\code{match.data} outputs matched data sets from
\code{matchit()}.
}

\usage{
  match.data(object, group="all", distance = "distance",
weights = "weights", subclass = "subclass")
}

\arguments{
  \item{object}{The output object from \code{matchit}. This is
    a required input.} 
  \item{group}{This argument specifies for which matched group the user
    wants to extract the data. Available options are \code{"all"} (all
    matched units), \code{"treat"} (matched units in the treatment
    group), and \code{"control"} (matched units in the control
    group). The default is \code{"all"}.}
  \item{distance}{This argument specifies the variable name used to
    store the distance measure. The default is \code{"distance"}.}
  \item{weights}{This argument specifies the variable name used to store
    the resulting weights from matching. The default is \code{"weights"}.}
  \item{subclass}{This argument specifies the variable name used to store the
    subclass indicator. The default is \code{"subclass"}.}
}

\value{
  Returns a subset of the original data set sent to
    \code{matchit()}, with just the matched units.  The data set also
    contains the additional variables \code{distance}, \code{weights},
    and \code{subclass}. The variable \code{distance} gives
    the estimated distance measure, and \code{weights} gives the
    weights for each unit, generated in the matching procedure.
    The variable \code{subclass} gives the subclass
    index for each unit (if applicable).
    See the \url{http://gking.harvard.edu/matchit/} for the complete
    documentation and type \code{demo(match.data)} at the R prompt to
    see a demonstration of the code.
}

\seealso{Please use \code{help.matchit} to access the matchit reference
  manual.  The complete document is available online at
  \url{http://gking.harvard.edu/matchit}.}

\author{
  Daniel Ho \email{daniel.ho@yale.edu}; Kosuke Imai
  \email{kimai@princeton.edu}; Gary King
  \email{king@harvard.edu}; Elizabeth Stuart \email{estuart@jhsph.edu}
}

\keyword{methods}
