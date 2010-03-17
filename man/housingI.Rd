\name{housingI}
\docType{data}
\alias{housingI}
\title{Boston Housing Interaction Data}
\description{
Median house price for 506 census tracts of Boston from the 1970 census.
The original data comprises 506 observations and 13 variables but has
been modified here to include all pairwise interactions of main effects
and to include B-spline basis functions of up to 6 degrees of freedom
for all original predictors.  In addition, all real valued variables
were mapped to dummy variables representing a factor with three levels
and all pairwise interactions of these dummy variables were added to the
design matrix.  In total, the modified data contains 506 observations
and 658 variables.  The outcome is the median house price \code{medv}.
}
\source{
Harrison D. and Rubinfeld D.L. (1978).
Hedonic prices and the demand for clean air.
\emph{J. Envir. Economics Management}, 5:81-102
}
\examples{data(housingI, package = "spikeslab")}
\keyword{datasets}
