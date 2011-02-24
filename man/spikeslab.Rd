\name{spikeslab}
\alias{spikeslab}
\title{Spike and Slab Regression}
\description{  
 Fits a rescaled spike and slab model using a continuous bimodal prior.
 A generalized elastic net estimator is used for variable selection and
 estimation.  Can be used for prediction and variable selection in low-
 and high-dimensional linear regression models.
}
\usage{
spikeslab(formula, data = NULL, x = NULL, y = NULL,
    n.iter1 = 500, n.iter2 = 500, mse = TRUE,
    bigp.smalln = FALSE, bigp.smalln.factor = 1, screen = (bigp.smalln),
    r.effects = NULL, max.var = 500, center = TRUE, intercept = TRUE,
    fast = TRUE, beta.blocks = 5, verbose = FALSE, ntree = 300,
    seed = NULL, ...)
}
\arguments{
    \item{formula}{A symbolic description of the model to be fit.}
    \item{data}{Data frame containing the data used in the formula.}
    \item{x}{x predictor matrix (can be used in place of formula and
      data frame call).}
    \item{y}{y response (can be used in place of formula and data frame call).}
    \item{n.iter1}{Number of burn-in Gibbs sampled values (i.e., discarded values).}
    \item{n.iter2}{Number of Gibbs sampled values, following burn-in.}
    \item{mse}{If TRUE, an external estimate for the
      overall variance is calculated using ridge regression or random
      forests (the latter is used when the degrees of freedom are low).
      Otherwise, the variance is included in the prior and estimated
      using Gibbs sampling.}
    \item{bigp.smalln}{Use if \code{p} >> \code{n}.}
    \item{bigp.smalln.factor}{Removes all variables except the top
      \code{n}times \code{bigp.smalln.factor} ones (used in filtering when 
      \code{p} >> \code{n}).}
    \item{screen}{If TRUE, variables are pre-filtered.}
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
\details{
  ---> General:
  
  The spike and slab method is described in detail in Ishwaran and Rao
  (2003, 2005a, 2005b and 2009).  For high-dimensional problems in which
  \code{p} >> \code{n}, where \code{p} is the number of variables and
  \code{n} is the sample size, use the option \code{bigp.smalln}=TRUE.
  Doing so implements a three-stage procedure:

    (1) Filtering step.  This removes all variables except the top
    \code{n} times \code{bigp.smalln.factor} ones.  Uses spike and slab
    regression with grouped regularization (complexity) parameters.
  
    (2) Model averaging step.  Refit the model using only those
    predictors from step 1.  Returns the posterior mean values from
    fitting a spike and slab model; referred to as the Bayesian model
    averaged (bma) estimate.
    
    (3) Variable selection step.  Select variables using the generalized
    elastic net (gnet).
    
  The filtering step is omitted when \code{bigp.smalln}=FALSE.
  Filtering can however be requested by setting \code{screen}=TRUE
  although users should be aware that this may degrade performance and
  should only be used when \code{p} is on the same order of \code{n}.

  Variables can be grouped using \code{r.effects}.  Grouping has the
  effect of forcing variables within a given group to share a common
  complexity (regularization) parameter.  To do so, define a list with
  each entry in the list made up of the variable names to be grouped.
  There is no limit to the number of groups.  Any variable that does
  not appear in the list will be assigned to a default group (the
  default group also has its own group-specific regularization
  parameter).  See Examples 1 and 3 below.

  
  ---> Miscellanea:

  By default, \code{fast}=TRUE when \code{bigp.smalln}=TRUE.  This
  invokes an ultra-fast filtering step.  Setting \code{fast}=FALSE
  invokes a more thorough filtering method that may slightly improve
  inferential results, but computational times will become very slow.
  The trade-off is unlikely to be justified.

  The formula and data-frame call should be avoided in high-dimensional
  problems and instead the x-predictor matrix and y response vector
  should be passed directly (see Example 3).  This avoids the huge
  overhead in parsing formula in R.

  By default, predictors are normalized to have mean 0 and variance 1.
  Pre-processing also involves centering y unless the user specifically
  requests that the intercept be excluded from the model.  Users can
  also over-ride centering predictors by setting \code{center}=FALSE.
  Use with extreme care.

  The \code{verbose} option sends output to the terminal showing the
  number of Gibbs iterations and the current complexity (regularization)
  parameter(s).
  
  Depends on the \code{randomForest} package for estimating the variance
  when \code{mse}=TRUE.  Note that \code{mse} is over-ridden and set to
  FALSE when \code{bigp.smalln}=TRUE.

  Depends on the \code{lars} package for the variable slection step.
}   
\value{
     An object of class \code{spikeslab} with the following components:
     \item{summary}{Summary object.}
     \item{verbose}{Verbose details (used for printing).}
     \item{terms}{Terms.}
     \item{sigma.hat}{Estimated variance.}
     \item{y}{Original y.}
     \item{xnew}{Centered, rescaled x-matrix.}
     \item{x}{Original x-matrix.}
     \item{y.center}{Centering for original y.}
     \item{x.center}{Centering for original x-matrix.}
     \item{x.scale}{Scaling for original x-matrix.}
     \item{names}{Variable names.}
     \item{bma}{bma coefficients in terms of xnew.}
     \item{bma.scale}{bma coefficients rescaled in terms of original x.}
     \item{gnet}{gnet coefficients in terms of xnew.}
     \item{gnet.scale}{gnet coefficients rescaled in terms of original x.}
     \item{gnet.path}{gnet path scaled in terms of the original x.}
     \item{gnet.obj}{gnet object (a lars-type object).}
     \item{gnet.obj.vars}{Variables (in order) used to calculate the gnet object.}
     \item{gnet.parms}{Generalized ridge regression parameters used to define the gnet.}
     \item{phat}{Estimated model dimension.}
     \item{complexity}{Complexity (regularization) parameter estimates.}
     \item{ridge}{List containing ridge values used to determine the bma.}
     \item{models}{List containing the models sampled.}
}
\author{
  Hemant Ishwaran (\email{hemant.ishwaran@gmail.com})

  J. Sunil Rao (\email{rao.jsunil@gmail.com})

  Udaya B. Kogalur (\email{kogalurshear@gmail.com})
}
\references{
  Breiman L. (2001). Random forests, \emph{Machine Learning}, 45:5-32.

  Efron B., Hastie T., Johnstone I. and Tibshirani R. (2004).
  Least angle regression (with discussion). \emph{Ann. Statist.},
  32:407-499.

  Ishwaran H. and Rao J.S. (2003).  Detecting differentially expressed
  genes in microarrays using Bayesian model selection.
  \emph{J. Amer. Stat. Assoc.}, 98:438-455.
  
  Ishwaran H. and Rao J.S. (2005a).  Spike and slab variable selection:
  frequentist and Bayesian strategies.  \emph{Ann. Statist.},
  33:730-773.

  Ishwaran H. and Rao J.S. (2005b).  Spike and slab gene selection for
  multigroup microarray data. \emph{J. Amer. Stat. Assoc.}, 100:764-780.

  Ishwaran H. and Rao J.S. (2010).  Generalized ridge regression:
  geometry and computational solutions when p is larger than n.

  Ishwaran H., Kogalur U.B. and Rao J.S. (2010). spikeslab: prediction
  and variable selection using spike and slab regression. \emph{R Journal},
  2(2), 68-73.

  Ishwaran H. and Rao J.S. (2011).  Mixing generalized ridge
  regressions.
  
  Zou H. and Hastie T. (2005).  Regularization and variable selection
  via the elastic net.  \emph{J. Royal Statist. Society B},
  67(2):301-320.
}
\seealso{
  \code{cv.spikeslab},
  \code{plot.spikeslab},
  \code{predict.spikeslab},
  \code{print.spikeslab},
  \code{sparsePC.spikeslab}.
}
\examples{

#------------------------------------------------------------
# Example 1:  diabetes data
#------------------------------------------------------------

# basic call
data(diabetesI, package = "spikeslab")
obj <- spikeslab(Y ~ . , diabetesI, verbose=TRUE)
print(obj)
plot(obj)

# grouping effect
# separate main effects and interactions into two groups
# use a group-specific regularization parameter for each group
xnames <- names(diabetesI[, -1])
r.eff <- vector("list", 2)
r.eff[[1]] <- xnames[c(1:10)]
r.eff[[2]] <- xnames[-c(1:10)]
obj2 <- spikeslab(Y ~ . , diabetesI, verbose=TRUE, r.effects=r.eff)
obj2
# extract the regularization parameters
print(apply(obj2$complexity, 2, summary))

\dontrun{
#------------------------------------------------------------
# Example 2: high-dimensional noise (diabetes data)
#------------------------------------------------------------

# add 2000 noise variables
data(diabetesI, package = "spikeslab")
diabetes.noise <- cbind(diabetesI,
      noise = matrix(rnorm(nrow(diabetesI) * 2000), nrow(diabetesI)))

# example of a big p, small n call
# don't use formula call; make call with x and y arguments
x <- diabetes.noise[, -1]
y <- diabetes.noise[, 1]
obj <- spikeslab(x=x, y=y, verbose=TRUE, bigp.smalln=TRUE, max.var=100)
obj

# same example ... but now group variables 
r.eff <- vector("list", 2)
r.eff[[1]] <- names(x)[c(1:100)]
r.eff[[2]] <- names(x)[-c(1:100)]
obj2 <- spikeslab(x=x, y=y, verbose=TRUE, bigp.smalln=TRUE,
                 r.effects=r.eff, max.var=100)
obj2

#------------------------------------------------------------
# Example 3: housing data with interactions
#------------------------------------------------------------

# another example of a big p, small n call
data(housingI, package = "spikeslab")
obj <- spikeslab(medv ~ ., housingI, verbose = TRUE,
           bigp.smalln = TRUE, max.var = 200)
print(obj)


}
}
\keyword{regression}

