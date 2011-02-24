####**********************************************************************
####**********************************************************************
####
####  SPIKE AND SLAB 1.1.3.1
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

plot.spikeslab <- function(x, plot.type = c("path", "cv"), breaks = FALSE, ...)
{

  ### Check that object is compatible
  if (!inherits(x, "spikeslab"))
     stop("This function only works for objects of class `spikeslab'")
  
  ###check whether object inherits mixing type
  ###make suitable alterations to merge mixing results
  if (sum(inherits(x, c("spikeslab", "mixing"), TRUE) == c(1, 2)) == 2) {
     x <- x$spikeslab.obj
  }

  ### determine the plot type
  plot.type <- match.arg(plot.type)

  ### plots
  if (plot.type == "path") {
    # plot the gnet solution path (calls plot.lars)
    x$gnet.obj$type <- "gnet"
    plot.lars(x$gnet.obj, breaks = breaks, ...)
  }
  else if (plot.type == "cv" & sum(inherits(x, c("spikeslab", "cv"), TRUE) == c(1, 2)) == 2) {
    # plot the cv curve
    cv <- x$cv
    cv.plot.path <- x$cv.path
    model.size <- unlist(x$gnet.path$model.size)
    K <- length(cv)
    p <- nrow(cv.plot.path) - 1
    cv.plot.mean <- apply(cv.plot.path, 1, mean, na.rm = TRUE)
    cv.plot.se <- apply(cv.plot.path, 1, SD)/sqrt(K)
    matplot(0:p, cv.plot.path, type = c("l", "n")[1 + 1 * (K > 20)], lty = 3, col = "gray", 
         xlim = range(c(0, model.size), na.rm = TRUE),
         ylim = range(c(cv.plot.path, cv.plot.mean + cv.plot.se, cv.plot.mean - cv.plot.se), na.rm = TRUE),
         xlab="Model Size", ylab="Cross-Validated MSE")
    lines(0:p, cv.plot.mean, lty = 1, lwd = 2, col = 4)
    error.bars(0:p, cv.plot.mean + cv.plot.se, cv.plot.mean - cv.plot.se, width = 0.0025, col = 2)
  }

}
