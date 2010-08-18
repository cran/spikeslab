\name{predict.spikeslab}
\alias{predict.spikeslab}
\title{Spike and Slab Prediction}
\description{
    Prediction on test data using spike and slab regression.
}
\usage{
\method{predict}{spikeslab}(object, newdata = NULL, ...)
}
\arguments{
  \item{object}{An object of class \code{spikeslab}.}
  \item{newdata}{Data frame or x-matrix containing test data (if
    omitted, the training data is used).}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Computes the predicted value using a test data set.
}    
\value{
  A vector of fitted values for the BMA and gnet and a matrix of fitted
  values for the gnet path.
}
\author{
  Hemant Ishwaran (\email{hemant.ishwaran@gmail.com})
  
  J. Sunil Rao (\email{rao.jsunil@gmail.com})
  
  Udaya B. Kogalur (\email{kogalurshear@gmail.com})
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
# Example 1: get the predictor for the training data
#------------------------------------------------------------
data(diabetesI, package = "spikeslab")
x <- diabetesI[, -1]
y <- diabetesI[, 1]
obj <- spikeslab(x = x, y = y)
#gnet predictor
yhat.gnet <- predict(obj)$yhat.gnet
#an equivalent call is...
yhat.gnet <- predict(obj, x = x)$yhat.gnet

#------------------------------------------------------------
# Example 2: ozone data with interactions
#------------------------------------------------------------

data(ozoneI, package = "spikeslab")
train.pt <- sample(1:nrow(ozoneI), nrow(ozoneI) * 0.80)
obj <- spikeslab(ozone ~ . , ozoneI[train.pt, ])
ytest <- ozoneI$ozone[-train.pt]
ss.pred <- predict(obj, ozoneI[-train.pt, ])
yhat.bma <- ss.pred$yhat.bma
yhat.gnet <- ss.pred$yhat.gnet
plot(ytest, yhat.bma, ylab = "yhat", pch = 16, col = 4)
points(ytest, yhat.gnet, pch = 16, col = 2)
abline(0, 1, lty = 2, col = 2)
legend("bottomright", legend = c("bma", "gnet"), col = c(4, 2), pch = 16)
}
}
\keyword{regression}
