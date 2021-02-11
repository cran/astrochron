\name{mwStats}
\alias{mwStats}
\title{'Dynamic window' moving average, median and variance of stratigraphic series}
\description{
'Dynamic window' moving average, median and variance of stratigraphic series. 
This routine adjusts the number of data points in the window so it has a constant duration in time or space, for use with unevenly sampled data.
}
\usage{
mwStats(dat,cols=NULL,win=NULL,conv=1,ends=F,CI=0,output=T,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Your data frame containing stratigraphic data; any number of columns (variables) are permitted, but the first column should be a location identifier (e.g., depth, height, time).}
  \item{cols}{A vector that identifies the variable column to be extracted (first column automatically extracted).}
  \item{win}{Moving window size in units of space or time.}
  \item{conv}{Convention for window placement: (1) center each window on a stratigraphic level in 'dat' (DEFAULT), (2) start with the smallest location datum in 'dat', (3) start with the largest location datum in 'dat'. For options 2 and 3, the center of the window will not necessarily coincide with a measured stratigraphic level in 'dat', but ends of the data set are better preserved. See options 'ends'.}
  \item{ends}{Assign average values to ends, by averaging data before first window, and averaging data after last window? (T or F; only applicable for conv=1)}
  \item{CI}{What confidence interval should be calculated for the average value (0-100 percent). If set to 0, the confidence interval calculation is skipped.}
  \item{output}{Output results? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}
\details{
If conv=1 is selected, the edges of the record are determined using a smaller window size. 
A constant value is assigned based on the observed values within the first and last 0.5*win 
of the record.
}

\value{
A data frame with five or six columns: 
Center of window, Average, Median, Variance, Number of points in window
If CI>0, the sixth column is the value used to determine the confidence interval (add and subtract it from the average) 
}

\examples{
# generate example series from ar1 noise, 5 kyr sampling interval
ex = ar1(npts=1001,dt=5)

# jitter sampling times
ex[1]=ex[1]+rnorm(1001,sd=1)
# sort
ex = ex[order(ex[,1],na.last=NA,decreasing=FALSE),]

# run mwStats
mwStats(ex,win=100)
}