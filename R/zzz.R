####**********************************************************************
####**********************************************************************
####
####  SPIKE AND SLAB 1.1.6
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
####  Written and Developed by:
####  ----------------------------------------------------------------
####  Hemant Ishwaran, Ph.D.
####  Director of Statistical Methodology
####  Professor, Division of Biostatistics
####  Clinical Research Building, Room 1058
####  1120 NW 14th Street
####  University of Miami, Miami FL 33136
####
####  email:  hemant.ishwaran@gmail.com
####  URL:    https://ishwaran.org/
####  ----------------------------------------------------------------
####  Udaya B. Kogalur, Ph.D.
####  Principal, Kogalur & Company, Inc.
####  5425 Nestleway Drive, Suite L
####  Clemmons, NC 27012
####
####  email:  ubk@kogalur.com
####  --------------------------------------------------------------
####
####**********************************************************************

.onAttach <- function(libname, pkgname) {
    ss.ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), 
                      fields="Version")
    packageStartupMessage(paste("\n",
                                pkgname,
                                ss.ver,
                                "\n",
                                "\n",
                                "Type spikeslab.news() to see new features, changes, and bug fixes.",
                                "\n",
                                "\n"))
}
