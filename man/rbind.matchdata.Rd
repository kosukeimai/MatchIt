\name{rbind.matchdata}
\alias{rbind.matchdata}
\alias{rbind.getmatches}

\title{
Append matched datasets together
}
\description{
These functions are \fun{rbind} methods for objects resulting from calls to \fun{match.data} and \fun{get_matches}. They function nearly identically to \code{rbind.data.frame()}; see Details for how they differ.
}
\usage{
\method{rbind}{matchdata}(\dots, deparse.level = 1)

\method{rbind}{getmatches}(\dots, deparse.level = 1)
}
\arguments{
  \item{\dots}{
Two or more \code{matchdata} or \code{getmatches} objects the output of calls to \fun{match.data} and \fun{get_matches}, respectively. Supplied objects must either be all \code{matchdata} objects or all \code{getmatches} objects.
}
  \item{deparse.level}{
Passed to \fun{rbind}.
}
}
\details{
\code{rbind()} appends two or more datasets row-wise. This can be useful when matching was performed separately on subsets of the original data and they are to be combined into a single dataset for effect estimation. Using the regular \code{data.frame} method for \code{rbind()} would pose a problem, however; the \code{subclass} variable would have repeated names across different datasets, even though units only belong to the subclasses in their respective datasets. \code{rbind.matchdata()} renames the subclasses so that the correct subclass membership is maintained.

The supplied matched datasets must be generated from the same original dataset, that is, having the same variables in it. The added components (e.g., weights, subclass) can be named differently in different datasets but will be changed to have the same name in the output.

\code{rbind.getmatches()} and \code{rbind.matchdata()} are identical.
}
\value{
An object of the same class as those supplied to it (i.e., a \code{matchdata} object if \code{matchdata} objects are supplied and a \code{getmatches} object if \code{getmatches} objects are supplied). \fun{rbind} is called on the objects after adjusting the variables so that the appropriate method will be dispatched corresponding to the class of the original data object.
}

\author{
Noah Greifer
}

\seealso{
\fun{match.data}, \fun{rbind}

See \code{vignettes("estimating-effects")} for details on using \code{rbind()} for effect estimation after subsetting the data.
}
\examples{
data("lalonde")

# Matching based on race subsets
m.out_b <- matchit(treat ~ age + educ + married +
                    nodegree + re74 + re75,
                  data = subset(lalonde, race == "black"))
md_b <- match.data(m.out_b)

m.out_h <- matchit(treat ~ age + educ + married +
                    nodegree + re74 + re75,
                  data = subset(lalonde, race == "hispan"))
md_h <- match.data(m.out_h)

m.out_w <- matchit(treat ~ age + educ + married +
                    nodegree + re74 + re75,
                  data = subset(lalonde, race == "white"))
md_w <- match.data(m.out_w)

#Bind the datasets together
md_all <- rbind(md_b, md_h, md_w)

#Subclass conflicts are avoided
levels(md_all$subclass)
}