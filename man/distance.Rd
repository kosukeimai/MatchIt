\name{distance}

\alias{distance}

\title{Distance function: estimating propensity scores}

\description{ The distance function calculates the distance measure to be used in the matching,
usually the propensity score. It is a sub-function of \code{matchit}.}

\details{This is a sub-function of the \code{matchit} command, which calculates the distance
measure used in the matching, usually the propensity score.  This function is called directly by
\code{matchit} and does not generally need to be called directly by users; these details are
included for advanced users. }


\usage{
distance <- distance(formula, model="logit", data, discard=0, reestimate=FALSE, counter=TRUE, ...)
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

	\item{model}{Method used to estimate the propensity score.  May be "logit" (default),
"probit", "nnet", "GAM", or "cart".}

	\item{discard}{Whether to discard units that fall outside some measure of support of the
distance score.  0 (default)=keep all units.  1=keep all units with common support.  2=discard
only control units outside the support of the distance measure of the treated units.  3=discard
only treated units outside the support of the distance measure of the control units.}

	\item{reestimate}{Specifies whether to reestimate the propensity score model after
discarding units (default=FALSE).}

	\item{counter}{Whether to display counter indicating the progress of the matching
(default=TRUE).}

	\item{...}{Additional arguments to be passed to \code{distance}, depending on the model
to be used.} }

\value{
  \item{in.sample}{Vector of length n showing whether each unit was eligible for matching due to
common support restrictions with \code{discard}.}
  \item{pscore}{Vector of estimated propensity scores.}
  \item{treat}{The treatment indicator from \code{data}.}
  \item{covariates}{Covariates used in the right-hand side of the assignment model.}
  \item{assign.model}{Output of the assignment model.}
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








