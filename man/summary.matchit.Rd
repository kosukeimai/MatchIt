\name{summary.matchit}

\alias{summary.matchit}
\alias{print.summary.matchit}

\title{Summarizing Quantities of Interest}

\description{ The \code{matchit} method for the generic \code{summary}
  command generates default output.
}

\usage{
\method{summary}{matchit}(object, verbose=F, sig=2, ...)
}

\arguments{ 
\item{object}{stored output from \code{matchit}.} 

\item{sig}{Optionally specifies the size of the t-statistic for a covariate to be labeled
"problematic" (default=2).}

\item{verbose}{Option to show summary statistics in \code{sum.all} and \code{sum.matched} for
all covariates, their squares, and two-way interactions when \code{verbose=TRUE} and only the
covariates themselves when \code{verbose=FALSE} (default).}

\item{\dots}{Additional parameters
passed to \code{summary.default}.} }

\value{ 
\item{call}{The original assignment model call.} 

\item{sum.all}{Data frame containing summary statistics on all observations in the treated and
control groups (means, standard deviations, t-statistics and bias statistics for all
covariates).}

\item{sum.matched}{Data frame containing summary statistics on the matched observations in the
treated and control groups (means, standard deviations, t-statistics, bias statistics, and
whether bias was reduced by the matching).}

\item{Problematic}{Problematic covariates that remain imbalanced in the matched samples.  
Variables listed when the absolute value of the t-statistic in the matched sample exceeds
\code{sig}.}

\item{xn}{Sample sizes in the full and matched samples.}

\item{q.table}{Array that contains the same information as \code{sum.all} and \code{sum.matched}
for each subclass.} 

\item{qn}{Sample sizes in the full and matched samples, by subclass.}

\item{match.matrix}{From the \code{matchit} output.} 

\item{psclass}{From the \code{matchit} output.} 
\item{in.sample}{From the \code{matchit} output.} }

\seealso{Please use \code{help.matchit} to access the matchit reference
  manual.  The complete document is available online at
  \url{http://gking.harvard.edu/matchit}.}

\author{
  Daniel Ho <\email{daniel.ho@yale.edu}>;  Kosuke Imai <\email{kimai@princeton.edu}>; Gary King
  <\email{king@harvard.edu}>; Elizabeth Stuart<\email{stuart@stat.harvard.edu}>
}

\keyword{methods}
