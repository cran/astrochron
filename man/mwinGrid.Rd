\name{mwinGrid}
\alias{mwinGrid}
\title{Determine 'dynamic moving window' for stratigraphic series, adjusting for changing sample density to maintain a window of constant duration; output on evenly spaced grid}
\description{
Determine start and end points for a moving window of fixed duration (e.g. 500 kiloyears). 
The dynamic window allows for adjustment of the number of data points in the window, 
so it has a constant duration in time or space. This version will output an evenly spaced
spatial/temporal grid.
}
\usage{
mwinGrid(dat,win,step,start=NULL,end=NULL,verbose=T)
}
\arguments{
  \item{dat}{Your data frame containing stratigraphic data; any number of columns (variables) are permitted, but the first column should be a location identifier (e.g., depth, height, time).}
  \item{win}{Moving window size in units of space or time.}
  \item{step}{Step size for moving window, in units of space or time.}
  \item{start}{Start moving window at what depth/height/time; by default will use first value}
  \item{end}{End moving window at what depth/height/time; by default will use last value}
  \item{verbose}{Verbose output? (T or F)}
}

\value{
A data frame containing: Starting index for window, Ending index for window, Location (center)
}

\details{
This algorithm is similar to function mwin, but instead of stepping forward one stratigraphic datum at a time, it generates an evenly spaced spatial/temporal grid.
}

\examples{
# generate some noise
ex1 <- ar1(npts=50,dt=1)

# jitter sampling times
ex1[1]=ex1[1]+rnorm(50,sd=0.25)
# sort 
ex1 = ex1[order(ex1[,1],na.last=NA,decreasing=FALSE),]

# run mwin
mwinGrid(ex1,win=10,step=2)
}