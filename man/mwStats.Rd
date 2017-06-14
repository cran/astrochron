\name{mwStats}
\alias{mwStats}
\title{'Dynamic window' moving average, median and standard deviation of stratigraphic series}
\description{
'Dynamic window' moving average, median and standard deviation of stratigraphic series. 
This routine adjusts the number of data points in the window so it has a constant duration in time or space, for use with unevenly sampled data.
}
\usage{
mwStats(dat,cols=NULL,win=NULL,conv=1,output=T,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Your data frame containing stratigraphic data; any number of columns (variables) are permitted, but the first column should be a location identifier (e.g., depth, height, time).}
  \item{cols}{A vector that identifies the variable column to be extracted (first column automatically extracted).}
  \item{win}{Moving window size in units of space or time.}
  \item{conv}{Convention for window placement: (1) center each window on a stratigraphic level in 'dat' (DEFAULT), (2) start with the smallest location datum in 'dat', (3) start with the largest location datum in 'dat'. For options 2 and 3, the center of the window will not necessarily coincide with a measured stratigraphic level in 'dat', but edges of the data set are better preserved.}
  \item{output}{Output results? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}
\value{
A data frame with four columns: Middle of window, Average, Median, Variance
}

\examples{
\donttest{
# generate example series from ar1 noise, 5 kyr sampling interval
ex = ar1(npts=1001,dt=5)

# jitter sampling times
ex[1]=ex[1]+rnorm(1001,sd=50)
# sort
ex = ex[order(ex[1],na.last=NA,decreasing=FALSE),]

# run mwStats
mwStats(ex,win=100)
}
}