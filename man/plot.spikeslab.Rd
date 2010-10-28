\name{plot.spikeslab}
\alias{plot.spikeslab}
\title{Plots for Spike and Slab Analysis}
\description{
    Plots either the gnet solution path or the cross-validated
    mean-squared-error (the latter only applies when cross-validation
    is used).
}
\usage{
\method{plot}{spikeslab}(x, plot.type = c("path", "cv"), breaks = FALSE, ...)
}
\arguments{
  \item{x}{An object of class \code{spikeslab}.}
  \item{plot.type}{Choosing "path" produces a plot of the gnet solution
       path. The choice "cv" produces the mean-squared error plot.  The
       latter applies only to objects from a \code{cv.spikeslab} call.}
  \item{breaks}{If TRUE, then vertical lines are drawn at each break
        point in the gnet solution path.}
  \item{...}{Further arguments passed to or from other methods.}
}
\author{
  Hemant Ishwaran (\email{hemant.ishwaran@gmail.com})
  
  J. Sunil Rao (\email{rao.jsunil@gmail.com})
  
  Udaya B. Kogalur (\email{kogalurshear@gmail.com})
}
\references{
  Efron B., Hastie T., Johnstone I., and Tibshirani R. (2004).
  Least angle regression (with discussion). \emph{Ann. Statist.},
  32:407-499.
   
  Ishwaran H. and Rao J.S. (2010).  Generalized ridge regression:
  geometry and computational solutions when p is larger than n.
}
\seealso{
  \code{spikeslab, cv.spikeslab}.
}
\examples{
\dontrun{
#------------------------------------------------------------
# Example 1: diabetes data
#------------------------------------------------------------

data(diabetesI, package = "spikeslab")
obj <- spikeslab(Y ~ . , diabetesI, verbose = TRUE)
plot(obj, plot.type = "path")
}
}
\keyword{plot}
