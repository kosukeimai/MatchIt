\name{match.data}

\alias{match.data}

\title{Output matched data sets}

\description{The code \code{match.data} creates output data sets from a \code{matchit}
matching algorithm.}

\usage{
match.data <- match.data(object, group="all")
}

\arguments{ 
\item{object}{(required).  Stored output from \code{matchit}.} 

\item{group}{Which units to output.  "all" (default) gives all matched units (treated and
control), "treat" gives just the matched treated units, and "control" gives just the matched
control units.} }

\value{ 
\item{Returns a subset of the original data set sent to \code{matchit}, with just the matched units.  The data set also contains the additional
variables \code{psclass}, \code{pscore}, and \code{psweights}.  The variable \code{psclass} gives the subclass index for each unit (if applicable).  
The variable \code{pscore} gives the propensity scores, and \code{psweights} gives the weights for each unit, generated in the matching procedure.  
See the \code{matchit} documentation for more details.} 
}

\seealso{Please use \code{help.matchit} to access the matchit reference
  manual.  The complete document is available online at
  \url{http://gking.harvard.edu/matchit}.}

\author{
  Daniel Ho <\email{deho@fas.harvard.edu}>;  Kosuke Imai <\email{kimai@princeton.edu}>; Gary King
  <\email{king@harvard.edu}>; Elizabeth Stuart<\email{stuart@stat.harvard.edu}>
}

\keyword{methods}
