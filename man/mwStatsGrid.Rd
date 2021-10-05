\name{mwStatsGrid}
\alias{mwStatsGrid}
\title{'Dynamic window' moving average, median and variance of stratigraphic series, using evenly spaced spatial/temporal grid}
\description{
'Dynamic window' moving average, median and variance of stratigraphic series. 
This routine adjusts the number of data points in the window so it has a constant duration in time or space, for use with unevenly sampled data.
The results are output on an evenly spaced spatial/temporal grid (this contrasts with mwStats).
}
\usage{
mwStatsGrid(dat,cols=NULL,win=NULL,step=NULL,start=NULL,end=NULL,output=T,norm=F,
            palette=6,ncolors=100,genplot=1,verbose=T)
}
\arguments{
  \item{dat}{Your data frame containing stratigraphic data; any number of columns (variables) are permitted, but the first column should be a location identifier (e.g., depth, height, time).}
  \item{cols}{A vector that identifies the variable column to be extracted (first column automatically extracted).}
  \item{win}{Moving window size, in units of space or time.}
  \item{step}{Moving window step size, in units of space or time.}
  \item{start}{Starting point for analysis, in units ofr space or time.}
  \item{end}{Ending point for analysis, in units of space or time.}
  \item{norm}{Normalize density estimates to maximum value? (T or F). If false, density estimates are normalized to unit area.}
  \item{output}{Output results? (T or F)}
  \item{palette}{What color palette would you like to use? (1) rainbow, (2) grayscale, (3) blue, (4) red, (5) blue-white-red, (6) viridis}
  \item{ncolors}{Number of colors to use in plot.}
  \item{genplot}{Generate summary plots? (0=none, 1=all time series, 2=kernel density estimates for each window, 3=kernel density estimates with median, 4=kernel density estimates with mean)}
  \item{verbose}{Verbose output? (T or F)}
}
\value{
A data frame with four columns: Center of window, Average, Median, Variance
}

\examples{
# generate example series from ar1 noise, 5 kyr sampling interval
ex = ar1(npts=1001,dt=5)

# jitter sampling times
ex[1]=ex[1]+rnorm(1001,sd=1)
# sort
ex = ex[order(ex[,1],na.last=NA,decreasing=FALSE),]

# run mwStats
mwStatsGrid(ex,win=100)
}