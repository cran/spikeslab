####**********************************************************************
####**********************************************************************
####
####  SPIKE AND SLAB 1.1.0
####
####  Copyright 2010, Cleveland Clinic Foundation
####
####  This program is free software; you can redistribute it and/or
####  modify it under the terms of the GNU General Public License
####  as published by the Free Software Foundation; either version 2
####  of the License, or (at your option) any later version.
####
####  This program is distributed in the hope that it will be useful,
####  but WITHOUT ANY WARRANTY; without even the implied warranty of
####  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
####  GNU General Public License for more details.
####
####  You should have received a copy of the GNU General Public
####  License along with this program; if not, write to the Free
####  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
####  Boston, MA  02110-1301, USA.
####
####  ----------------------------------------------------------------
####  Project Partially Funded By:
####    --------------------------------------------------------------
####    National Science Foundation, Grants DMS-0705037, DMS-0405675 and DMS-0405072
####
####    Hemant Ishwaran, Ph.D.
####    Dept of Quantitative Health Sciences/Wb4
####    Cleveland Clinic Foundation
####    9500 Euclid Avenue
####    Cleveland, OH 44195
####
####    email:  hemant.ishwaran@gmail.com
####    phone:  216-444-9932
####    URL:    www.bio.ri.ccf.org/Resume/Pages/Ishwaran/ishwaran.html
####
####
####	J. Sunil Rao, Ph.D.
####    Deparment of Biostatistics
####    University of Miami
####
####    email: rao.jsunil@gmail.com
####
####    --------------------------------------------------------------
####    Case Western Reserve University/Cleveland Clinic  
####    CTSA Grant:  XX1 RR000000, National Center for
####    Research Resources (NCRR), NIH
####
####  ----------------------------------------------------------------
####  Written by:
####    --------------------------------------------------------------
####    Hemant Ishwaran, Ph.D.
####    Dept of Quantitative Health Sciences/Wb4
####    Cleveland Clinic Foundation
####    9500 Euclid Avenue
####    Cleveland, OH 44195
####
####    email:  hemant.ishwaran@gmail.com
####    phone:  216-444-9932
####    URL:    www.bio.ri.ccf.org/Resume/Pages/Ishwaran/ishwaran.html
####
####  ----------------------------------------------------------------
####  Maintained by:
####    Udaya B. Kogalur, Ph.D.
####    Dept of Quantitative Health Sciences/Wb4
####    Cleveland Clinic Foundation
####    
####    Kogalur Shear Corporation
####    5425 Nestleway Drive, Suite L1
####    Clemmons, NC 27012
####
####    email:  kogalurshear@gmail.com
####    phone:  919-824-9825
####    URL:    www.kogalur-shear.com
####    --------------------------------------------------------------
####
####**********************************************************************
####**********************************************************************

cv.spikeslab <- function(
 x=NULL,              #x matrix
 y=NULL,              #y response
 K=10,                #K-fold
 parallel=FALSE,      #parallel processing via snow
 plot.it=TRUE,        #plot cv
 n.iter1=500,         #no. burn-in samples
 n.iter2=500,         #no. Gibbs sampled values (following burn-in)
 mse=TRUE,            #mse estimate (TRUE --> ridge/forest estimate)
 bigp.smalln=FALSE,   #used for p>>n 
 bigp.smalln.factor=1,#p>>n factor (relative to n) used in filtering variables
 screen=(bigp.smalln),#filter variables?
 r.effects=NULL,      #used for grouping variables
 max.var=500,         #max no. vars allowed in final model
 center=TRUE,         #center data (FALSE --> used for array data: !!CAUTION!!)
 intercept=TRUE,      #include an intercept?
 fast=TRUE,           #update beta in blocks (for bigp.small n, controls screening)
 beta.blocks=5,       #no. of beta blocks in beta Gibbs update (fast=TRUE)
 verbose=TRUE,        #verbose details
 ntree=300,           #number RF trees
 seed=NULL,           #seed
  ...)
{

#tidy up x
x <- as.matrix(x)
if (length(unique(colnames(x))) != ncol(x)) {
    colnames(x) <- paste("x.", 1:ncol(x), sep = "")
}
  
# define the folds
# last fold is the full data and corresponds to the "primary object"
all.folds <-  split(sample(1:nrow(x)), rep(1:K, length = nrow(x)))
K <- length(all.folds)
all.folds[[K+1]] <- nrow(x) + 1

#core cv function
eval.fold <- function(k, ...) {
  if (verbose) {
    if (k <= K) {
      cat("\t K-fold:", k, "\n")
    }
    else {
      cat("\t final analysis (full-data)\n")
    }
  }
  omit <- all.folds[[k]]
  obj <- spikeslab(x = as.matrix(x[-omit,, drop = FALSE]), y = y[-omit],
                     n.iter1 = n.iter1, n.iter2 = n.iter2,
                     mse = mse, bigp.smalln = bigp.smalln,
                     bigp.smalln.factor = bigp.smalln.factor,
                     screen = screen, r.effects = r.effects,
                     max.var = max.var, center = center, intercept = intercept,
                     fast = fast, beta.blocks = beta.blocks, verbose = FALSE,
                     ntree = ntree, seed = seed)
  if (k <= K) {
    #test-set prediction
    pred.obj <- predict(obj, as.matrix(x[omit,, drop = FALSE])) 
    yhat.k <- pred.obj$yhat.gnet
    yhat.path.k <- pred.obj$yhat.gnet.path
    gnet.path.k <- obj$gnet.path$path
    #lars should only return steps 0, ..., p; yet there seems to be ties
    #apply an ad-hoc beta-breaker here
    model.size.k <- apply(gnet.path.k, 1,
             function(sbeta){sum(abs(sbeta) > .Machine$double.eps, na.rm = TRUE)})
    beta.break.k <- which(!duplicated(model.size.k))
    if (length(beta.break.k) > 0) {
      gnet.path.k <- as.matrix(gnet.path.k[beta.break.k, ])
      yhat.path.k <- as.matrix(yhat.path.k[, beta.break.k])
      #test-set error
      cv.k <- mean((y[omit] - yhat.k)^2, na.rm = TRUE)
      cv.path.k <- colMeans((y[omit] - yhat.path.k)^2, na.rm = TRUE)
      gnet.k <- gnet.path.k[which(cv.path.k == min(cv.path.k))[1], ]
    }
    else {
      cv.path.k <- cv.k <- mean((y[omit] - mean(y[-omit]))^2, na.rm = TRUE)
      gnet.k <- rep(0, length(obj$names))
    }
    return(list(model.size.k = model.size.k,
              cv.k = cv.k,
              cv.path.k = cv.path.k,
              gnet.k = gnet.k,
              names = obj$names))
  }
  else {
    return(list(obj = obj))
  }
}

## Determine if parallel processing is to be done.
if (parallel) {
 if (!require(snow, quietly = TRUE)) {
     warning("package 'snow' not found, i.e., parallel processing was not performed")
     eval.fold.obj <- lapply(1:(K+1), eval.fold, ...)
  }   
  else {
    ## The default number of threads is two (2).
    if (parallel == TRUE) {
      parallel <- 2
    }
    cl <- makeSOCKcluster(rep("localhost", as.integer(parallel)))

    test <-  list(x, y, n.iter1, n.iter2, mse, bigp.smalln, bigp.smalln.factor,
                  screen, r.effects, max.var, center, intercept,
                  fast, beta.blocks, verbose, ntree, seed, all.folds)

    clusterEvalQ(cl, library(spikeslab))
    
    eval.fold.obj <- clusterApplyLB(cl, 1:(K+1), eval.fold, test)

    stopCluster(cl)
  }
}
else {
  eval.fold.obj <- lapply(1:(K+1), eval.fold, ...)
}


#extract the primary obj
#redefine eval.fold.obj
primary.obj <- eval.fold.obj[[K+1]]$obj
eval.fold.obj[[K+1]] <- NULL
varnames <- primary.obj$names
p <- length(varnames)

#parse the eval.fold.obj
cv <- model.size <- rep(NA, K)
cv.path <- list(length = K)
gnet.path <- stability <- matrix(0, K, p)
cv.plot.path <-  matrix(NA, p + 1, K)
for (k in 1:K) {
  #extract cv, cv.path, gnet, gnet.path, stability, model size
  cv[k] <- eval.fold.obj[[k]]$cv.k
  cv.path[[k]] <- eval.fold.obj[[k]]$cv.path.k
  gnet.k <- eval.fold.obj[[k]]$gnet.k
  gnet.path[k, is.element(varnames, eval.fold.obj[[k]]$names)] <- gnet.k
  stability[k, is.element(varnames, eval.fold.obj[[k]]$names)] <- 1 * (abs(gnet.k) >  .Machine$double.eps)
  model.size[k] <- max(eval.fold.obj[[k]]$model.size.k)
  #convert mse into matrix format more conducive for plotting/printing
  cv.plot.path[1:length(cv.path[[k]]), k] <- cv.path[[k]]
} 
cv.plot.mean <- apply(cv.plot.path, 1, mean, na.rm = TRUE)
cv.plot.se <- apply(cv.plot.path, 1, SD)/sqrt(K)

## plot it
if (plot.it) {
  matplot(0:p, cv.plot.path, type = c("l", "n")[1 + 1 * (K > 20)], lty = 3, col = "gray", 
         xlim = range(c(0, model.size), na.rm = TRUE),
         ylim = range(c(cv.plot.path, cv.plot.mean + cv.plot.se,
            cv.plot.mean - cv.plot.se), na.rm = TRUE),
         xlab="Model Size", ylab="Cross-Validated MSE")
  lines(0:p, cv.plot.mean, lty = 1, lwd = 2, col = 4)
  error.bars(0:p, cv.plot.mean + cv.plot.se, cv.plot.mean - cv.plot.se, width = 0.0025, col = 2)
}

# stability analysis; make it pretty for the return
tally.stability <- cbind(primary.obj$gnet,
                         apply(gnet.path, 2, mean, na.rm = TRUE) * primary.obj$x.scale,
                         primary.obj$gnet.scale,
                         apply(gnet.path, 2, mean, na.rm = TRUE),
                         100 * apply(stability, 2, mean, na.rm = TRUE))
colnames(tally.stability) <- c("gnet", "gnet.cv", "gnet.scale", "gnet.scale.cv", "stability")
rownames(tally.stability) <- varnames
tally.stability <- tally.stability[order(tally.stability[, 5], abs(tally.stability[, 1]),
                          decreasing = TRUE),, drop = FALSE]

# pretty details for terminal output
get.model.size <- function(mn, se) {
  ms.upper <- min(which(mn == min(mn, na.rm = TRUE)))
  mn.upper <- mn[ms.upper]
  ms.range <- which((mn + se >= mn.upper) & (mn - se <= mn.upper))
  ms.lower <- min(ms.range)
  ms.all <- unique(c(ms.lower, ms.upper))
  if (length(ms.all) == 1) return(ms.all)
  paste("[", ms.lower, ",", ms.upper, "]", sep = "")
}
if (verbose){
cat("-------------------------------------------------------------------","\n")
cat("Big p small n                 :",bigp.smalln,"\n")
cat("Screen variables              :",screen,"\n")
cat("Fast processing               :",fast,"\n")
cat("Sample size                   :",nrow(x),"\n")
cat("No. predictors                :",ncol(x),"\n")
cat("No. burn-in values            :",n.iter1,"\n")
cat("No. sampled values            :",n.iter2,"\n")
cat("K-fold                        :",K,"\n")
cat("CV mean-squared error         :", paste(round(mean(cv, na.rm = TRUE), 3) , "+/-",
                                         round(sd(cv, na.rm = TRUE)/sqrt(K), 3), "\n"))
cat("Model size                    :",get.model.size(cv.plot.mean, cv.plot.se), "\n")
cat("-------------------------------------------------------------------","\n")
}

#return the goodies
object <- list(cv = cv, cv.path = cv.plot.path, stability = tally.stability,
               gnet.path = primary.obj$gnet.path, gnet.obj = primary.obj$gnet.obj)
invisible(object)

}
