\name{plot.matchit}

\alias{plot.matchit}

\title{Graphing Quantities of Interest}

\description{ The \code{matchit} method for the generic \code{plot}
  command generates default plots.}

\usage{
\method{plot}{matchit}(x, ...)
}

\arguments{
\item{x}{stored output from \code{\link{matchit}}.  }
\item{\dots}{Additional parameters passed to \code{plot.default}.
  Because \code{plot.matchit} primarily produces diagnostic plots, many
  of these parameters are hard-coded for convenience and
  presentation. }  
}

\value{ 
\item{\code{plot.matchit} allows the user to check the distribution of all covariates
in the assignment model, squares, and interactions, within all subclasses.  The graphs present
density estimate graphs of the propensity score of treated and control units in the full and
matched samples, jitter plots of the propensity score for treated and control units, density
estimate graphs of any covariates, and density estimate graphs of any covariates by subclass.}}

\seealso{Please use \code{help.matchit} to access the matchit reference
  manual.  The complete document is available online at
  \url{http://gking.harvard.edu/matchit}.
}
  
\author{
  Daniel Ho <\email{deho@fas.harvard.edu}>;  Kosuke Imai <\email{kimai@princeton.edu}>; Gary King
  <\email{king@harvard.edu}>; Elizabeth Stuart<\email{stuart@stat.harvard.edu}>
}

\keyword{hplot}
