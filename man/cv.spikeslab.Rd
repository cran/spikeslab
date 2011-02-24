\name{cv.spikeslab}
\alias{cv.spikeslab}
\title{K-fold Cross-Validation for Spike and Slab Regression}
\description{
  Computes the K-fold cross-validated mean squared prediction error for
  the generalized elastic net from spike and slab regression.  Returns a
  stability index for each variable.
}
\usage{
cv.spikeslab(x = NULL, y = NULL, K = 10, parallel = FALSE, 
    plot.it = TRUE, n.iter1 = 500, n.iter2 = 500, mse = TRUE,
    bigp.smalln = FALSE, bigp.smalln.factor = 1, screen = (bigp.smalln),
    r.effects = NULL, max.var = 500, center = TRUE, intercept = TRUE,
    fast = TRUE, beta.blocks = 5, verbose = TRUE, save.all = TRUE,
    ntree = 300, seed = NULL, ...)
}
\arguments{
    \item{x}{x-predictor matrix.}
    \item{y}{y-response values.}
    \item{K}{Number of folds.}
    \item{parallel}{Indicates whether computations should be performed
      using parallel processing.  Parallel processing is implemented via
      the package \code{snow}.  This package should be loaded prior to
      calling \code{cv.spikeslab}.  When not set to FALSE, the variable
      should indicate the number of parallel threads the user desires to
      create.  Currently, parallel execution is limited to socket
      clusters on the local machine on which the R session is being
      initiated.  Thus, a machine with 4 CPU's will allow one thread to
      be dedicated to each CPU by setting \code{parallel} to 4.  Setting
      \code{parallel} to a number greater than the number of CPU's on
      the machine is not recommended as this will not increase
      performance.  The default number of parallel threads is 2.}
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
    \item{save.all}{If TRUE, spikeslab object for each fold is saved and returned.}
    \item{ntree}{Number of trees used by random forests (applies only when \code{mse} is TRUE).}
    \item{seed}{Seed for random number generator.  Must be a negative
        integer.}
    \item{...}{Further arguments passed to or from other methods.}
}
\value{
  Invisibly returns a list with components:
    \item{spikeslab.obj}{Spike and slab object from the full data.}
    \item{cv.spikeslab.obj}{List containing spike and slab objects from each fold.
      Can be NULL.}
    \item{cv.fold}{List containing the cv splits.}
    \item{cv}{Mean-squared error for each fold for the gnet.}
    \item{cv.path}{A matrix of mean-squared errors for the gnet
       solution path. Rows correspond to model sizes, columns are the folds.}
    \item{stability}{Matrix containing stability for each variable defined as the
       percentage of times a variable is identified over the K-folds.
       Also includes bma and gnet coefficient values and their cv-fold-averaged
       values.}
    \item{bma}{bma coefficients from the full data in terms of the standardized x.}
    \item{bma.scale}{bma coefficients from the full data, scaled in terms of the original x.}
    \item{gnet}{cv-optimized gnet in terms of the standardized x.}
    \item{gnet.scale}{cv-optimized gnet in terms of the original x.}
    \item{gnet.model}{List of models selected by gnet over the K-folds.}
    \item{gnet.path}{gnet path from the full data, scaled in terms of the original x.}
    \item{gnet.obj}{gnet object from fitting the full data (a lars-type object).}
    \item{gnet.obj.vars}{Variables (in order) used to calculate the gnet object.}
    \item{verbose}{Verbose details (used for printing).}
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

  Ishwaran H. and Rao J.S. (2010).  Generalized ridge regression:
  geometry and computational solutions when p is larger than n.

  Ishwaran H. and Rao J.S. (2011).  Mixing generalized ridge
  regressions.
}
\seealso{
  \code{sparsePC.spikeslab},
  \code{plot.spikeslab},
  \code{predict.spikeslab},
  \code{print.spikeslab}.
}
\examples{
\dontrun{
#------------------------------------------------------------
# Example 1: 10-fold validation using parallel processing
#------------------------------------------------------------

data(ozoneI, package = "spikeslab")
y <- ozoneI[,  1]
x <- ozoneI[, -1]
cv.obj <- cv.spikeslab(x = x, y = y, parallel = 4)
plot(cv.obj, plot.type = "cv")
plot(cv.obj, plot.type = "path")

#------------------------------------------------------------
# Example 2: 10-fold validation using parallel processing
# (high dimensional diabetes data)
#------------------------------------------------------------

# add 2000 noise variables
data(diabetesI, package = "spikeslab")
diabetes.noise <- cbind(diabetesI,
      noise = matrix(rnorm(nrow(diabetesI) * 2000), nrow(diabetesI)))
x <- diabetes.noise[, -1]
y <- diabetes.noise[, 1]

cv.obj <- cv.spikeslab(x = x, y = y, bigp.smalln=TRUE, parallel = 4)
plot(cv.obj)
}
}
\keyword{regression}

