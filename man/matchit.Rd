\name{matchit}

\alias{matchit}
\alias{MatchIt}
\alias{Matchit}

\title{MatchIt: Matching Software for Causal Inference}

\description{ \code{matchit} is the main command of the package
\emph{MatchIt}, which enables parametric models for causal inference to
work better by selecting well-matched subsets of the original treated
and control groups.  MatchIt implements the suggestions of Ho, Imai,
King, and Stuart (2004) for improving parametric statistical models by
preprocessing data with nonparametric matching methods.  MatchIt
implements a wide range of sophisticated matching methods, making it
possible to greatly reduce the dependence of causal inferences on
hard-to-justify, but commonly made, statistical modeling assumptions.
The software also easily fits into existing research practices since,
after preprocessing with MatchIt, researchers can use whatever
parametric model they would have used without MatchIt, but produce
inferences with substantially more robustness and less sensitivity to
modeling assumptions.  Matched data sets created by MatchIt can be
entered easily in Zelig (\url{http://gking.harvard.edu/zelig}) for
subsequent parametric analyses. Full documentation is available online
at \url{http://gking.harvard.edu/matchit}, and help for specific
commands is available through \code{help.matchit}.}

\details{The matching is done using the \code{matchit(treat ~ X, ...)} 
command, where \code{treat} is the vector of treatment assignments and
\code{X} are the covariates to be used in the matching.  There are a
number of matching options, detailed below.  The full syntax is
\code{matchit(formula, data=NULL, discard=0, exact=FALSE, replace=FALSE,
ratio=1, model="logit", reestimate=FALSE, nearest=TRUE, m.order=2,
caliper=0, calclosest=FALSE, mahvars=NULL, subclass=0, sub.by="treat",
counter=TRUE, full=FALSE, full.options=list(), \dots)} A summary of the
results can be seen graphically using \code{plot(matchitobject)}, or
numerically using \code{summary(matchitobject)}.
\code{print(matchitobject)} also prints out the output.  }

\usage{matchit(formula, data, method = "nearest", distance = "logit",
               distance.options = list(), discard = "none",
               reestimate = FALSE, ...)
}


\arguments{
  \item{formula}{This argument takes the usual syntax of R formula,
    \code{treat ~ x1 + x2}, where \code{treat} is a binary treatment
    indicator and \code{x1} and \code{x2} are the pre-treatment
    covariates. Both the treatment indicator and pre-treatment covariates
    must be contained in the same data frame, which is specified as
    \code{data} (see below).  All of the usual R syntax for formula
    works. For example, \code{x1:x2} represents the first order
    interaction term between \code{x1} and \code{x2}, and \code{I(x1^2)}
represents the square term of \code{x1}. See \code{help(formula)}
    for details.}
  \item{data}{This argument specifies the data frame containing the
    variables called in \code{formula}.}
  \item{method}{This argument specifies a matching method. Currently,
    \code{"exact"} (exact matching), \code{"full"} (full matching),
    \code{"genetic"} (genetic matching), \code{"nearest"} (nearest
    neighbor matching), \code{"optimal"} (optimal matching), and
    \code{"subclass"} (subclassification) are available. The default is
    \code{"nearest"}. Note that within each of these matching methods,
    \emph{MatchIt} offers a variety of options.}
  \item{distance}{This argument specifies the method used to estimate the
    distance measure. The default is logistic regression,
    \code{"logit"}. A variety of other methods are available.}
  \item{distance.options}{ This optional argument specifies the optional
    arguments that are passed to the model for estimating the distance
    measure. The input to this argument should be a list.}
  \item{discard}{This argument specifies whether to discard units that
    fall outside some measure of support of the distance score before
    matching, and not allow them to be used at all in the matching
    procedure.  Note that discarding units may change the quantity of
    interest being estimated.  The options are: \code{"none"} (default), which discards no units before matching,
    \code{"both"}, which discards all units (treated and control) that are
      outside the support of the distance measure,
    \code{"control"}, which discards only control units outside the
      support of the distance measure of the treated units, and
    \code{"treat"}, which discards only treated units outside the support
      of the distance measure of the control units.}
  \item{reestimate}{This argument specifies whether the model for
    distance measure should be re-estimated after units are
    discarded. The input must be a logical value. The default is
    \code{FALSE}.}
  \item{...}{Additional arguments to be passed to a variety of matching
    methods.}
}

\value{
  \item{call}{The original \code{matchit} call.}
  \item{formula}{The formula used to specify the model for
    estimating the distance measure.}
  \item{model}{The output of the model used to estimate
    the distance measure.  \code{summary(m.out$model)} will give the
    summary of the model where \code{m.out} is the output object from
    \code{matchit}.}
  \item{match.matrix}{An \eqn{n_1} by \code{ratio} matrix
    where the row names, which can be obtained through
    \code{row.names(match.matrix)}, represent the names of the
    treatment units, which come from the data frame specified in
    \code{data}.  Each column stores the name(s) of the control unit(s) matched
    to the treatment unit of that row. For example, when the
    \code{ratio} input for nearest neighbor or optimal matching is
    specified as 3, the three columns of 
    \code{match.matrix} represent the three control units matched to
    one treatment unit).    
    \code{NA} indicates that the treatment unit was not matched.}   
  \item{discarded}{A vector of length $n$ that displays
    whether the units were ineligible for matching due to common
    support restrictions.  It equals \code{TRUE} if unit \eqn{i} was
    discarded, and it is set to \code{FALSE} otherwise.}
  \item{distance}{A vector of length \eqn{n} with the estimated
    distance measure for each unit.}
  \item{weights}{A vector of length \eqn{n} that provides the
    weights assigned to each unit in the matching process.  Unmatched
    units have weights equal to \code{0}. Matched treated units have
    weight \code{1}.  Each matched control unit has weight proportional
    to the number of treatment units to which it was matched, and the sum of
    the control weights is equal to the number of uniquely matched
    control units.}
  \item{subclass}{The subclass index in an ordinal
    scale from 1 to the total number of subclasses as specified in
    \code{subclass} (or the total number of subclasses from full or
    exact matching).  Unmatched units have \code{NA}.}
  \item{q.cut}{The subclass cut-points that classify the
    distance measure.}
  \item{treat}{The treatment indicator from
    \code{data} (the left-hand side of \code{formula}).}
  \item{X}{The covariates used for estimating the
    distance measure (the right-hand side of \code{formula}).}
  \item{nn}{A basic summary table of matched data (e.g., the number of matched units)}
}

\seealso{Please use \code{help.matchit} to access the matchit reference
  manual.  The complete document is available online at
  \url{http://gking.harvard.edu/matchit}.  
}

\references{Daniel Ho, Kosuke Imai, Gary King, and Elizabeth Stuart
  (2007).  Matching as Nonparametric Preprocessing for Reducing Model Dependence in Parametric
  Causal Inference.  Political Analysis 15(3): 199-236.
  \url{http://gking.harvard.edu/files/abs/matchp-abs.shtml} 
}

\author{
  Daniel Ho \email{daniel.ho@yale.edu};  Kosuke Imai \email{kimai@princeton.edu}; Gary King
  \email{king@harvard.edu}; Elizabeth
  Stuart\email{estuart@jhsph.edu} }

\keyword{environment}








