\name{print.spikeslab}
\alias{print.spikeslab}
\title{Print Summary Output of Analysis}
\description{
    Print summary output from spike and slab analysis.
    Note that this is the default print method for the package.
}
\usage{
\method{print}{spikeslab}(x, ...)
}
\arguments{
    \item{x}{An object of class \code{spikeslab}.}
    \item{...}{Further arguments passed to or from other methods.}
}
\author{
  Hemant Ishwaran (\email{hemant.ishwaran@gmail.com})
  
  J. Sunil Rao (\email{rao.jsunil@gmail.com})
  
  Udaya B. Kogalur (\email{ubk@kogalur.com})
}
\references{
  Ishwaran H. and Rao J.S. (2003).  Detecting differentially expressed
  genes in microarrays using Bayesian model selection.
  \emph{J. Amer. Stat. Assoc.}, 98:438-455.
  
  Ishwaran H. and Rao J.S. (2005a).  Spike and slab variable selection: frequentist and Bayesian
  strategies.  \emph{Ann. Statist.}, 33:730-773.

  Ishwaran H. and Rao J.S. (2005b).  Spike and slab gene selection for
  multigroup microarray data. \emph{J. Amer. Stat. Assoc.}, 100:764-780.

  Ishwaran H. and Rao J.S. (2009).  Generalized ridge regression:
  geometry and computational solutions when p is larger than n.
}
\seealso{
  \code{spikeslab}.
}
\examples{
\dontrun{
#------------------------------------------------------------
# Example 1: diabetes data
#------------------------------------------------------------

data(diabetesI, package = "spikeslab")
obj <- spikeslab(Y ~ . , diabetesI, verbose = TRUE)
print(obj)
}
}
\keyword{print}
