\name{diagnose}

\alias{diagnose}

\title{Diagnostics for matching procedure}

\description{\code{diagnose} is a sub-function of \code{matchit} which calculates summary
statistics for the matching, including the weight for each unit.}


\details{This is a sub-function of \code{matchit} which calculates summary statistics for the
matching procedure in \code{matchit}, calculating t-statistics for covariates as well as
weights for each unit.  This function is called directly by \code{matchit} and does not
generally need to be called directly by users; these details are included for advanced users. }

\usage{ diagnose <- diagnose(formula, match.matrix, pscore, in.sample, data, exact=FALSE,
		mahvars=NULL, subclass=0, psclass=NULL, nearest=TRUE, q.cut=NULL, counter=TRUE)
}


\arguments{
  \item{formula}{(required).  Takes the form of \code{T ~ X1 + X2}, where \code{T} is a binary
treatment indicator and \code{X1} and \code{X2} are the pre-treatment covariates, and \code{T},
\code{X1}, and \code{X2} are contained in the same data frame.  The \code{+} symbol means
"inclusion" not "addition." You may also include interaction terms in the form if
\code{I(X1*X2)} or squared terms in the form of \code{I(X1^2)}.}

  \item{match.matrix}{(required). n1 by ratio data frame where the rows correspond to treated
units and the columns store the names of the control units matched to each treated unit.  NA
indicates that treated unit was not matched. Generally created in \code{nearest}.}

  \item{pscore}{(required). Vector of estimated propensity scores.  Generally calculated in
\code{distance}.}

  \item{in.sample}{(required). Vector of length n showing whether each unit was eligible for
matching due to common support restrictions with \code{discard}.  Generally calculated in
\code{distance}.}

  \item{data}{(required).  Data frame containing the variables called in the \code{formula}.  
The dataframe should not include variables with the names \code{psclass}, \code{psweights}, or
\code{pscore}, as these are expressly reserved in the output dataframe for MatchIt.}

	\item{exact}{"FALSE" (default)=no exact matching.  "TRUE"=exact matching on all
variables in \code{formula}.  A vector of variable names (that are in \code{data} to indicate
separate variables on which to exact match, in combination with matching on the propensity
score.}

	\item{mahvars}{Variables on which to perform Mahalanobis matching within each caliper
(default=NULL).  Should be entered as a vector of names of variables in \code{data}.}

	\item{subclass}{Either a scaler specifying the number of subclasses (default=0) or a
vector of probabilities to create quantiles based on \code{sub.by}.}

  \item{psclass}{Subclass index in an ordinal scale from 1 to the number of subclasses.  
Unmatched units have subclass=0.  Generally computed in \code{subclassify}.}

  \item{nearest}{Whether to perform nearest-neighbor matching (default=TRUE).  }
  \item{q.cut}{Subclass cut points.  Generally calculated in \code{subclassify}.}

   \item{counter}{Whether to display counter indicating the progress of the matching
(default=TRUE).} }

\value{
  \item{psweights}{Vector of length n giving the weight assigned to each unit in the matching
process.  Each weight is proportional to the number of times that unit was matched.} }

\seealso{Please use \code{help.matchit} to access the matchit reference
  manual.  The complete document is available online at
  \url{http://gking.harvard.edu/matchit}.  
}

\author{
  Daniel Ho <\email{deho@fas.harvard.edu}>;  Kosuke Imai <\email{kimai@princeton.edu}>; Gary King
  <\email{king@harvard.edu}>; Elizabeth Stuart<\email{stuart@stat.harvard.edu}>
}

\keyword{methods}








