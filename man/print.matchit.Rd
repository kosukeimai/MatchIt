\name{print.matchit}

\alias{print.matchit}

\title{Printing Quantities of Interest}

\description{ The \code{matchit} method for the generic \code{print}
  command generates default output.
}

\usage{
\method{print}{matchit}(x,  digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
\item{x}{stored output from \code{matchit}. } 
\item{digits}{Optionally specify the significant digits to be displayed.}
\item{\dots}{Additional parameters passed to \code{print.default}.}  
}

\value{ 
\item{Returns the original assignment model call, summary statistics of the propensity
score for the full and matched samples, summary statistics of the propensity score for each
subclass (if applicable), and sample sizes for the full and matched samples (full and by
subclass).} }

\seealso{Please use \code{help.matchit} to access the matchit reference
  manual.  The complete document is available online at
  \url{http://gking.harvard.edu/matchit}.}

\author{
  Daniel Ho <\email{deho@fas.harvard.edu}>;  Kosuke Imai <\email{kimai@princeton.edu}>; Gary King
  <\email{king@harvard.edu}>; Elizabeth Stuart<\email{stuart@stat.harvard.edu}>
}

\keyword{methods}
