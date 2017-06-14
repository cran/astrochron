\name{mwin}
\alias{mwin}
\title{Determine 'dynamic moving window' for stratigraphic series}
\description{
Determine start and end points for a moving window of fixed duration (e.g. 500 kiloyears). 
The dynamic window allows for adjustment of the number of data points in the window, 
so it has a constant duration in time or space.
}
\usage{
mwin(dat,win,conv=1,verbose=T)
}
\arguments{
  \item{dat}{Your data frame containing stratigraphic data; any number of columns (variables) are permitted, but the first column should be a location identifier (e.g., depth, height, time).}
  \item{win}{Moving window size in units of space or time.}
  \item{conv}{Convention for window placement: (1) center each window on a stratigraphic level in 'dat' (DEFAULT), (2) start with the smallest location datum in 'dat', (3) start with the largest location datum in 'dat'. For options 2 and 3, the center of the window will not necessarily coincide with a measured stratigraphic level in 'dat', but edges of the data set are better preserved.}
  \item{verbose}{Verbose output? (T or F)}
}

\value{
A data frame containing: Starting index for window, Ending index for window, Location (average), Location (center), Location (midpoint)
}

\details{
This algorithm steps forward one stratigraphic datum at a time. The output consists of:

Average = this is the average of the depth/time values in the given window.

Center = this is the center of the 'win' size window.

Midpoint = this is midpoint between first and last observation in the window 
(for unevenly sampled data this is typically less than than the size of 'win').
}

\examples{
\donttest{
# generate some noise
ex1 <- ar1(npts=50,dt=1)

# jitter sampling times
ex1[1]=ex1[1]+rnorm(50,sd=3)
# sort 
ex1 = ex1[order(ex1[1],na.last=NA,decreasing=FALSE),]

# run mwin
mwin(ex1,win=10)
}
}