\name{cv.spikeslab}
\alias{cv.spikeslab}
\title{K-fold Cross-Validation for Spike and Slab Regression}
\description{
  Computes the K-fold cross-validated mean squared prediction error for
  the generalized elastic net from spike and slab regression.  Returns a
  stability index for each variable.
}
\usage{
cv.spikeslab(x = NULL, y = NULL, K = 10, plot.it = TRUE,
    n.iter1 = 500, n.iter2 = 500, mse = TRUE,
    bigp.smalln = FALSE, bigp.smalln.factor = 1, screen = (bigp.smalln),
    r.effects = NULL, max.var = 500, center = TRUE, intercept = TRUE,
    fast = TRUE, beta.blocks = 5, verbose = TRUE, ntree = 300,
    seed = NULL, ...)
}
\arguments{
    \item{x}{x-predictor matrix.}
    \item{y}{y-response values.}
    \item{K}{Number of folds.}
    \item{plot.it}{If TRUE, plots the mean prediction error and its standard error.}
    \item{n.iter1}{Number of burn-in Gibbs sampled values (i.e., discarded values).}
    \item{n.iter2}{Number of Gibbs sampled values, following burn-in.}
    \item{mse}{If TRUE, an external estimate for the overall variance is calculated.} 
    \item{bigp.smalln}{Use if \code{p} >> \code{n}.}
    \item{bigp.smalln.factor}{Top \code{n} times this value of variables
      to be kept in the filtering step (used when \code{p} >> \code{n}).}
    \item{screen}{If TRUE, variables are first pre-filtered.}
    \item{r.effects}{List used for grouping variables (see details below).}
    \item{max.var}{Maximum number of variables allowed in the final model.}
    \item{center}{If TRUE, variables are centered by their
      means. Default is TRUE and should only be adjusted in extreme examples.}
    \item{intercept}{If TRUE, an intercept is included in the model,
      otherwise no intercept is included.  Default is TRUE.}
    \item{fast}{If TRUE, use blocked Gibbs sampling to accelerate the algorithm.}
    \item{beta.blocks}{Update beta using this number of blocks (\code{fast}
      must be TRUE).}
    \item{verbose}{If TRUE, verbose output is sent to the terminal.}
    \item{ntree}{Number of trees used by random forests (applies only when \code{mse} is TRUE).}
    \item{seed}{Seed for random number generator.  Must be a negative
        integer.}
    \item{...}{Further arguments passed to or from other methods.}
}
\value{
     Invisibly returns a list with components:
     \item{cv}{K-dimensional vector of mean-squared errors for the gnet.}
     \item{cv.path}{A matrix of mean-squared errors for the gnet
       solution path.  Rows correspond to model sizes, columns are the folds.}
     \item{stability}{Stability for each variable defined as the
       percentage of times a variable is identified over the K-folds.}
     \item{gnet.path}{gnet path from the full data, scaled in terms of the original x.}
     \item{gnet.obj}{gnet object from fitting the full data (a lars-type object).}
}
\author{
  Hemant Ishwaran (\email{hemant.ishwaran@gmail.com})

  J. Sunil Rao (\email{rao.jsunil@gmail.com})

  Udaya B. Kogalur (\email{kogalurshear@gmail.com})
}
\references{
  Ishwaran H. and Rao J.S. (2005a).  Spike and slab variable selection:
  frequentist and Bayesian strategies.  \emph{Ann. Statist.},
  33:730-773.

  Ishwaran H. and Rao J.S. (2009).  Generalized ridge regression:
  geometry and computational solutions when p is larger than n.
}
\seealso{
  \code{sparsePC.spikeslab},
  \code{predict.spikeslab},
  \code{print.spikeslab}.
}
\examples{
\dontrun{
#------------------------------------------------------------
# Example:  10-fold validation
#------------------------------------------------------------

data(diabetesI, package = "spikeslab")
y <- diabetesI[,  1]
x <- diabetesI[, -1]
cv.obj <- cv.spikeslab(x = x, y = y)
print(head(cv.obj$stability, 25))
}
}
\keyword{regression}

