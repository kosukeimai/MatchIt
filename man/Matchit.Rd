\name{Matchit}

\alias{Matchit}

\title{Matchit: Matching Software for Causal Inference}

\description{

\emph{Matchit} enables parametric models for causal inference to work better by selecting
well-matched subsets of the original treated and control groups.  MatchIt implements the
suggestions of Ho, Imai, King, and Stuart (2004) for improving parametric statistical models by
preprocessing data with nonparametric matching methods.  MatchIt implements a wide range of
sophisticated matching methods, making it possible to greatly reduce the dependence of causal
inferences on hard-to-justify, but commonly made, statistical modeling assumptions.  The
software also easily fits into existing research practices since, after preprocessing with
MatchIt, researchers can use whatever parametric model they would have used without MatchIt, but
produce inferences with substantially more robustness and less sensitivity to modeling
assumptions.  Matched data sets created by MatchIt can be entered easily in Zelig
(\url{http://gking.harvard.edu/zelig}) for subsequent parametric analyses. Full documentation is
available online at \url{http://gking.harvard.edu/matchit}, and help for specific commands is
available through \code{help.matchit}.}

\details{
The matching is done using the \code{matchit(treat ~ X, ...)} command, where \code{treat} is the
vector of treatment assignments and \code{X} are the covariates to be used in the matching.  
There are a number of matching options, detailed below.  The full syntax is
\code{matchit(formula, data=NULL, discard=0, exact=FALSE, replace=FALSE, ratio=1, model="logit",
reestimate=FALSE, nearest=TRUE, m.order=2, caliper=0, calclosest=FALSE, mahvars=NULL,
subclass=0, subtype=1, interact=TRUE, \dots)} A summary of the results can be seen graphically
using \code{plot(matchitobject)}, or numerically using \code{summary(matchitobject)}.  
\code{print(matchitobject)} also prints out the output.  }

\usage{
matchit <- matchit(formula, data, model="logit", discard=0, reestimate=FALSE, nearest=TRUE,
                 replace=FALSE, m.order=2, ratio=1, caliper=0, calclosest=FALSE,
                 subclass=0, sub.by="treat", mahvars=NULL, exact=FALSE, counter=TRUE, ...)
}


\arguments{ 
\item{formula}{(required).  Takes the form of \code{T ~ X1 + X2}, where \code{T} is
a binary treatment indicator and \code{X1} and \code{X2} are the pre-treatment covariates, and
\code{T}, \code{X1}, and \code{X2} are contained in the same data frame.  The \code{+} symbol
means "inclusion" not "addition." You may also include interaction terms in the form if
\code{I(X1*X2)} or squared terms in the form of \code{I(X1^2)}.}

\item{data}{(required).  Data frame containing the variables called in the \code{formula}.  The
dataframe should not include variables with the names \code{psclass}, \code{psweights}, or
\code{pscore}, as these are expressly reserved in the output dataframe for MatchIt.}
	
\item{model}{Method used to estimate the propensity score.  May be "logit" (default), "probit",
"nnet", "GAM", or "cart".}
	
\item{discard}{Whether to discard units that fall outside some measure of support of the
distance score.  0 (default)=keep all units.  1=keep all units with common support.  2=discard
only control units outside the support of the distance measure of the treated units.  3=discard
only treated units outside the support of the distance measure of the control units.}
	
\item{reestimate}{Specifies whether to reestimate the propensity score model after discarding
units (default=FALSE).}
	
\item{nearest}{Whether to perform nearest-neighbor matching (default=TRUE).  }
\item{replace}{Whether to match with replacement (default=FALSE). } \item{m.order}{Order in
which to match treated units with control units.  2 (default)=from high to low, 3=from low to
high, 4=random order.} 

\item{ratio}{The number of control units to be matched to each treated
unit (default=1).} 

\item{caliper}{Standard deviations of the propensity score within which to
draw control units (default=0).} 

\item{calclosest}{If \code{caliper!=0}, whether to take the
nearest available match if no matches are available within \code{caliper} (default=FALSE).}
	
\item{subclass}{Either a scaler specifying the number of subclasses (default=0) or a vector of
probabilities to create quantiles based on \code{sub.by}.}
	
\item{sub.by}{If \code{subclass!=0}, by what criteria to subclassify.  "treat" (default) =by the
number of treated units, "control"=by the number of control units, "all"=by the total number of
units.}
	
\item{mahvars}{Variables on which to perform Mahalanobis matching within each caliper
(default=NULL).  Should be entered as a vector of names of variables in \code{data}.}
	
\item{exact}{"FALSE" (default)=no exact matching.  "TRUE"=exact matching on all variables in
\code{formula}.  A vector of variable names (that are in \code{data} to indicate separate
variables on which to exact match, in combination with matching on the propensity score.}
	
\item{counter}{Whether to display counter indicating the progress of the matching
(default=TRUE).}
	
\item{...}{Additional arguments to be passed to \code{matchit}, depending on the model to be
used.} 
}

\value{
  \item{call}{The original \code{matchit} call.}
  \item{formula}{Formula used to specify the propensity score.}
  \item{match.matrix}{n1 by ratio data frame where the rows correspond to treated units and the
columns store the names of the control units matched to each treated unit.  NA indicates that
treated unit was not matched.}
 
 \item{in.sample}{Vector of length n showing whether each unit was eligible for matching due to
common support restrictions with \code{discard}.}

  \item{matched}{Vector of length n showing whether each unit was matched.}

  \item{psweights}{Vector of length n giving the weight assigned to each unit in the matching
process.  Each weight is proportional to the number of times that unit was matched.}

  \item{psclass}{Subclass index in an ordinal scale from 1 to the number of subclasses.  
Unmatched units have subclass=0.}
  \item{q.cut}{Subclass cut points.}
  \item{assign.model}{Output of the assignment model.}
  \item{data}{The original data set, with \code{psclass}, \code{psweights}, and \code{pscore}
(propensity scores) added.}
  \item{treat}{The treatment indicator from \code{data}.}
  \item{covariates}{Covariates used in the right-hand side of the assignment model.}
}

\seealso{Please use \code{help.matchit} to access the matchit reference
  manual.  The complete document is available online at
  \url{http://gking.harvard.edu/matchit}.  
}

\references{Daniel Ho, Kosuke Imai, Gary King, and Elizabeth Stuart (2004)
}

\author{
  Daniel Ho <\email{deho@fas.harvard.edu}>;  Kosuke Imai <\email{kimai@princeton.edu}>; Gary King
  <\email{king@harvard.edu}>; Elizabeth Stuart<\email{stuart@stat.harvard.edu}>
}

\keyword{environment}








