\name{sparsePC.spikeslab}
\alias{sparsePC.spikeslab}
\alias{sparsePC}
\title{Multiclass Prediction using Spike and Slab Regression}
\description{
 Variable selection for the multiclass gene prediction problem.
}
\usage{
sparsePC.spikeslab(x = NULL, y = NULL, n.rep = 10,
  n.iter1 = 150, n.iter2 = 100, n.prcmp = 5, max.genes = 100,
  ntree = 1000, nodesize = 1, verbose = TRUE, ...)
}
\arguments{
    \item{x}{x matrix of gene expressions.}
    \item{y}{Class labels.}
    \item{n.rep}{Number of Monte Carlo replicates.} 
    \item{n.iter1}{Number of burn-in Gibbs sampled values (i.e., discarded values).}
    \item{n.iter2}{Number of Gibbs sampled values, following burn-in.}
    \item{n.prcmp}{Number of principal components.}
    \item{max.genes}{Maximum number of genes in final model.}
    \item{ntree}{Number of trees used by random forests.}
    \item{nodesize}{Nodesize of trees.}    
    \item{verbose}{If TRUE, verbose output is sent to the terminal.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Multiclass prediction using a hybrid combination of spike and slab
  linear regression and random forest multiclass prediction (Ishwaran
  and Rao, 2009).  A pseudo y-vector of response values is calculated
  using each of the top \code{n.prcmp} principal components of the
  x-gene expression matrix.  The generalized elastic net obtained from
  using spike and slab regression is used to select genes; one
  regression fit is used for each of the pseduo y-response vectors.  The
  final combined set of genes are passed to random forests and used to
  construct a multiclass forest predictor.  This procedure is repeated
  \code{n.rep} times with each Monte Carlo replicate based on balanced
  cross-validation with 2/3rds of the data used for training and 1/3rd
  used for testing.

  ---> Miscellanea:

  Test set error is only computed when \code{n.rep} is larger than 1.
  If \code{n.rep}=1 the full data is used without any cross-validation.
}
\value{
  Invisibly, the final set of selected genes as well as the complete set
  of genes selected over the \code{n.rep} Monte Carlo replications.  The
  random forest classifier is also returned.

  The misclassification error rate, error rate for each class, and other
  summary information are output to the terminal.
}
\author{
  Hemant Ishwaran (\email{hemant.ishwaran@gmail.com})
  
  J. Sunil Rao (\email{rao.jsunil@gmail.com})
  
  Udaya B. Kogalur (\email{ubk@kogalur.com})
}
\references{
  Ishwaran H. and Rao J.S. (2009).  Generalized ridge regression:
  geometry and computational solutions when p is larger than n.
}
\seealso{
  \code{spikeslab}.
}
\examples{
\dontrun{
#------------------------------------------------------------
# Example 1:  leukemia data
#------------------------------------------------------------

data(leukemia, package = "spikeslab")
sparsePC.out <- sparsePC(x = leukemia[, -1], y = leukemia[, 1], n.rep = 3)
rf.obj <- sparsePC.out$rf.obj
varImpPlot(rf.obj)
}
}
\keyword{regression}
