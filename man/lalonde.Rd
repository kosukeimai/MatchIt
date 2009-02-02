\name{lalonde}
\docType{data}
\alias{lalonde}
\title{Data from National Supported Work Demonstration and PSID, as analyzed by Dehejia and Wahba (1999).}

\description{ This is a subsample of the data from the treated group in the National Supported
Work Demonstration (NSW) and the comparison sample from the Current Population Survey (CPS).  
This data was previously analyzed extensively by Lalonde (1986) and Dehejia and Wahba (1999).  
The full dataset is available at \url{http://www.columbia.edu/~rd247/nswdata.html}. }

\usage{data(lalonde)}

\format{ A data frame with 313 observations (185 treated, 429 control).  There are 10 variables
measured for each individual.  "treat" is the treatment assignment (1=treated, 0=control).  
"age" is age in years.  "educ" is education in number of years of schooling.  "black" is an
indicator for African-American (1=African-American, 0=not).  "hispan" is an indicator for being
of Hispanic origin (1=Hispanic, 0=not).  "married" is an indicator for married (1=married, 0=not
married).  "nodegree" is an indicator for whether the individual has a high school degree (1=no
degree, 0=degree). "re74" is income in 1974, in U.S. dollars.  "re75" is income in 1975, in U.S.
dollars.  "re78" is income in 1978, in U.S. dollars. }

\references{
	Lalonde, R. (1986). Evaluating the econometric evaluations of training programs 
	with experimental data. American Economic Review 76: 604-620.

	Dehejia, R.H. and Wahba, S. (1999).  Causal Effects in Nonexperimental Studies: 
	Re-Evaluating the Evaluation  of Training Programs.  Journal of the American 
	Statistical Association 94: 1053-1062.
}
\source{\url{http://www.columbia.edu/~rd247/nswdata.html}}
\keyword{datasets}  
