\name{exactmatch}

\alias{exactmatch}

\title{Exact matching}

\description{\code{exactmatch} is a sub-function of \code{matchit} that does exact matching
on the set of covariates.}


\details{This sub-function is called automatically by the \code{matchit} command when
\code{exact=TRUE} is chosen, meaning exact matching on all covariates in the assignment model.  
It creates subclasses defined by the values of the covariates, where each subclass must contain
some treated and some control units.  Units with no exact matches in the other treatment group
are discarded.  This function is called directly by \code{matchit} and does not generally need
to be called directly by users; these details are included for advanced users. }

\usage{
exactmatch <- exactmatch(formula, data, counter=TRUE)
}


\arguments{
  \item{formula}{(required).  Takes the form of \code{T ~ X1 + X2}, where \code{T} is a binary
treatment indicator and \code{X1} and \code{X2} are the pre-treatment covariates, and \code{T},
\code{X1}, and \code{X2} are contained in the same data frame.  Exact matching is done on the
variables on the right-hand side of \code{formula}.}

  \item{data}{(required).  Data frame containing the variables called in the \code{formula}.  
The dataframe should not include variables with the names \code{psclass}, \code{psweights}, or
\code{pscore}, as these are expressly reserved in the output dataframe for MatchIt.}

	\item{counter}{Whether to display counter indicating the progress of the matching
(default=TRUE).} }

\value{
  \item{psclass}{Subclass index in an ordinal scale from 1 to the number of subclasses.  
Unmatched units have subclass=0.  With exact matching, these subclasses are defined by the
categories of the variables in \code{formula}.} }

\seealso{Please use \code{help.matchit} to access the matchit reference
  manual.  The complete document is available online at
  \url{http://gking.harvard.edu/matchit}.  
}

\author{
  Daniel Ho <\email{deho@fas.harvard.edu}>;  Kosuke Imai <\email{kimai@princeton.edu}>; Gary King
  <\email{king@harvard.edu}>; Elizabeth Stuart<\email{stuart@stat.harvard.edu}>
}

\keyword{methods}








