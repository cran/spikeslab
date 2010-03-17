\name{ozoneI}
\docType{data}
\alias{ozoneI}
\title{Ozone Interaction Data}
\description{
The data consists of 366 readings of maximum daily ozone measured in the
Los Angeles basin.  After removing missing values, the original data has
been expanded to include all pairwise interactions, as well as B-spline
basis functions (6 degrees of freedom), for each of the original 12
variables (9 meteorlogical variables and 3 variables recording date of
measurement: month, day of the month, and day of week).  In total, the
modified data has 203 observations and 134 variables.  The outcome is
\code{ozone}.
}
\source{
Breiman L. and Friedman J.H. (1985).
Estimating optimal transformations for multiple regression and correlation.
\emph{J. Amer. Stat. Assoc.}, 80:580-598.
}
\examples{data(ozoneI, package = "spikeslab")}
\keyword{datasets}

