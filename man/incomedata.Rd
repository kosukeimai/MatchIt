\name{incomedata}
\docType{data}
\alias{incomedata}
\title{Simulated data to illustrate matching procedure for Matchit}

\description{
	This data set is a fake data set to be used to illustrate matching methods in the
Matchit library.  It has data on 200 individuals, 100 treated and 100 control, with treatment
assignemt, four covariates (income in 1974, income in 1975, male, and education), and one
outcome variable (income in 1978) observed. } 

\usage{data(incomedata)} 

\format{
	A data frame with 200 observations (rows) and six variables (columns).  The first
column, "treat", is the treatment assignment (1=treatment, 0=control).  The second column,
"income74", is income in 1974 (in U.S. dollars).  The third column, "income75", is income in
1975, also in U.S. dollars.  The fourth column, "male", is an indicator for male (1=male,
0=female). The fifth column, "educ", is the education level in years of schooling.  The sixth
column is the outcome variable, income in 1978 ("income78").  } 

\source{Simulated data.}
\keyword{datasets}
