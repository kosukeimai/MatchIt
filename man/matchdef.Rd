\name{matchdef}

\alias{matchdef}

\title{Nearest neighbor matching algorithm}

\description{ This is a sub-function of \code{matchit} which performs the nearest neighbor
matching; \code{matchdef} is called automatically by \code{matchit}. }


\details{The matching is done using nearest neighbor matching on the propensity score, estimated
in the sub-function \ref{distance}.  This function is called directly by \code{matchit} and
does not generally need to be called directly by users; these details are included for advanced
users. }

\usage{
matchdef <- matchdef(formula, in.sample, pscore, nearest=TRUE, replace=FALSE, 
		m.order=2, ratio=1, caliper=0, calclosest=FALSE, mahvars=NULL, exact=FALSE, 
		data=NULL, counter=TRUE)
}


\arguments{
  \item{formula}{(required).  Takes the form of \code{T ~ X1 + X2}, where \code{T} is a binary
treatment indicator and \code{X1} and \code{X2} are the pre-treatment covariates, and \code{T},
\code{X1}, and \code{X2} are contained in the same data frame.  The \code{+} symbol means
"inclusion" not "addition." You may also include interaction terms in the form if
\code{I(X1*X2)} or squared terms in the form of \code{I(X1^2)}.}

  \item{data}{Data frame containing the variables in \code{formula}.}
  \item{in.sample}{Vector of length n showing whether each unit was eligible for matching due to
common support restrictions with \code{discard}. Computed in \code{distance}.}

  \item{pscore}{Vector of propensity scores, calculated in \code{distance}.}

	\item{nearest}{Whether to perform nearest-neighbor matching (default=TRUE).  }
	\item{replace}{Whether to match with replacement (default=FALSE). }
	\item{m.order}{Order in which to match treated units with control units.  1=optimal (requires ``optmatch" package), 
2 (default)=from high to low, 3=from low to high, 4=random order.}
	\item{ratio}{The number of control units to be matched to each treated unit (default=1).}

	\item{caliper}{Standard deviations of the propensity score within which to draw control
units (default=0).}

	\item{calclosest}{If \code{caliper!=0}, whether to take the nearest available match if
no matches are available within \code{caliper} (default=FALSE).}

	\item{mahvars}{Variables on which to perform Mahalanobis matching within each caliper
(default=NULL).  Should be entered as a vector of names of variables in \code{data}.}

	\item{exact}{"FALSE" (default)=no exact matching.  "TRUE"=exact matching on all
variables in \code{formula}.  A vector of variable names (that are in \code{data}) to indicate
separate variables on which to exact match, in combination with matching on the propensity
score.}

	\item{counter}{Whether to display counter indicating the progress of the matching
(default=TRUE).} }

\value{
  \item{match.matrix}{n1 by ratio data frame where the rows correspond to treated units and the
columns store the names of the control units matched to each treated unit.  NA indicates that
treated unit was not matched.}

  \item{in.sample}{Vector of length n showing whether each unit was eligible for matching due to
common support restrictions with \code{discard}.} }

\seealso{Please use \code{help.matchit} to access the matchit reference
  manual.  The complete document is available online at
  \url{http://gking.harvard.edu/matchit}.  
}

\author{
  Daniel Ho <\email{daniel.ho@yale.edu}>;  Kosuke Imai <\email{kimai@princeton.edu}>; Gary King
  <\email{king@harvard.edu}>; Elizabeth Stuart<\email{stuart@stat.harvard.edu}>
}

\keyword{internal}








