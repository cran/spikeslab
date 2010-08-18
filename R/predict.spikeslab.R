####**********************************************************************
####**********************************************************************
####
####  SPIKE AND SLAB 1.1.1
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

predict.spikeslab <- function(object, newdata = NULL, ...)
{

 ### --------------------------------------------------------------
 ### preliminary checks

  #check that object is compatible
  if (!inherits(object,"spikeslab"))
     stop("This function only works for objects of class `spikeslab'")
  #check for newdata
  if (is.null(newdata)) {
    newdata <- object$x
    colnames(newdata) <- object$names
  }
  #special allowance for newdata with one observation
  newdata <- rbind(newdata)
  #check for colnames
  if (is.null(colnames(newdata))) {
    if (ncol(newdata) == length(object$names)) {
       colnames(newdata) <- object$names
    } else {
       stop("Number of columns of testing data does not match training data")
    }
  }
  
 ### --------------------------------------------------------------
 ### ensure coherence between test and training data

 #restrict predictors in test data to those in the training data
 newdata <- as.data.frame(rbind(newdata))
 if (is.null(object$terms)) {
   old.xnames <- object$names
 }
 else {
   old.xnames <- attr(object$terms, "term.labels")
 }
 new.xnames <- names(newdata)
 if (sum(!is.element(old.xnames, new.xnames)) > 0)
   stop("predictors in test data don't match training data")
 newdata <- newdata[, match(old.xnames, new.xnames), drop = FALSE]
 
 #build the test design matrix
 if (is.null(object$terms)) {
   x.new <- newdata
 }
 else {
   Terms <- delete.response(object$terms)
   old.na.action <- options()$na.action
   na.keep <- function(x){x}
   options(na.action = na.keep)
   x.new <- model.matrix(Terms, newdata)
   options(na.action=old.na.action)
 }

 #ensure that test data matches training data one more time
 #needed because of factors
 cov.names <- unlist(dimnames(x.new)[2])
 x.new <- x.new[, match(object$names, cov.names), drop = FALSE]

 #remove NA's 
 x.new <- as.matrix(rbind(x.new))
 has.NA <- which(apply(x.new,1,function(x){any(is.na(x))}))
 if (length(has.NA) > 0){
   x.new <- as.matrix(x.new[-has.NA,, drop = FALSE])
 }
 rownames(x.new) <- colnames(x.new) <- NULL
  
 ### --------------------------------------------------------------
 ### extract x/y centering information
 ### get coefficient values
 ### correct for division by zero
 ### NA coefficient estimates ---> 0
  
 y.center <- object$y.center
 x.center <- object$x.center
 bma.scale <- object$bma.scale
 bma.scale[is.na(bma.scale)] <- 0
 gnet.scale <- object$gnet.scale
 gnet.scale[is.na(gnet.scale)] <- 0
 gnet.path <- object$gnet.path$path
 gnet.path[is.na(gnet.path)] <- 0

 ### --------------------------------------------------------------
 ### predicted values

 yhat.bma <- c(y.center - sum(bma.scale * x.center) + x.new %*% bma.scale)
 yhat.gnet <- c(y.center - sum(gnet.scale * x.center) + x.new %*% gnet.scale)
 yhat.gnet.path <- apply(gnet.path, 1, function(sbeta) {
                         c(y.center - sum(sbeta * x.center) + x.new %*% sbeta)})
   
 ### --------------------------------------------------------------
 #### return the goodies

 return(list(yhat.bma = yhat.bma, yhat.gnet = yhat.gnet, yhat.gnet.path = yhat.gnet.path))

}
