\name{subclassify}

\alias{subclassify}

\title{Subclassify units}

\description{This is a sub-function of \code{matchit} which generates subclasses of units
given a propensity score.}

\details{This function creates subclasses of units based on their covariate values, usually
defined by the propensity score.  This function is called directly by \code{matchit} and does
not generally need to be called directly by users; these details are included for advanced
users. }

\usage{
subclassify <- subclassify(formula, data, in.sample, pscore, nearest=TRUE, match.matrix, 
	      	subclass=0, sub.by="treat", counter=TRUE, full=FALSE, full.options=list())
}


\arguments{
  \item{formula}{(required).  Takes the form of \code{T ~ X1 + X2}, where \code{T} is a binary
treatment indicator and \code{X1} and \code{X2} are the pre-treatment covariates, and \code{T},
\code{X1}, and \code{X2} are contained in the same data frame.  The \code{+} symbol means
"inclusion" not "addition." You may also include interaction terms in the form if
\code{I(X1*X2)} or squared terms in the form of \code{I(X1^2)}.}

  \item{data}{(required).  Data frame containing the variables called in the \code{formula}.  
The dataframe should not include variables with the names \code{psclass}, \code{psweights}, or
\code{pscore}, as these are expressly reserved in the output dataframe for MatchIt.}

  \item{in.sample}{(required).  Vector of length n showing whether each unit was eligible for
matching due to common support restrictions with \code{discard}.  Generally created by
\code{distance}.}

  \item{pscore}{(required).  Vector of estimated propensity scores.  Generally created by
\code{distance}.}

	\item{nearest}{Whether to perform nearest-neighbor matching (default=TRUE).  }

  \item{match.matrix}{(required) n1 by ratio data frame where the rows correspond to treated
units and the columns store the names of the control units matched to each treated unit.  NA
indicates that treated unit was not matched.  Generally created by \code{matchdef}.}

	\item{subclass}{Either a scaler specifying the number of subclasses (default=0) or a
vector of probabilities to create quantiles based on \code{sub.by}.}

	\item{sub.by}{If \code{subclass!=0}, by what criteria to subclassify.  "treat" (default)
=by the number of treated units, "control"=by the number of control units, "all"=by the total
number of units.}

	\item{counter}{Whether to display counter indicating the progress of the matching
(default=TRUE).} 

   \item{full}{Whether to do full matching (default=FALSE).  Requires ``optmatch" package.}
   \item{full.options}{Additional options for full matching.}


}

\value{
  \item{psclass}{Subclass index in an ordinal scale from 1 to the number of subclasses.  
Unmatched units have subclass=0.}
  \item{q.cut}{Subclass cut points.}
}

\seealso{Please use \code{help.matchit} to access the matchit reference
  manual.  The complete document is available online at
  \url{http://gking.harvard.edu/matchit}.  
}

\author{
  Daniel Ho <\email{deho@fas.harvard.edu}>;  Kosuke Imai <\email{kimai@princeton.edu}>; Gary King
  <\email{king@harvard.edu}>; Elizabeth Stuart<\email{stuart@stat.harvard.edu}>
}

\keyword{internal}








