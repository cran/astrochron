\name{mwCor}
\alias{mwCor}
\title{Calculate moving window correlation coefficient for two stratigraphic series, using a 'dynamic window'}
\description{
Calculate moving window correlation coefficient for two stratigraphic series, using a 'dynamic window'. 
This routine adjusts the number of data points in the window so it has a constant duration in time or space, for use with unevenly sampled data.
}
\usage{
mwCor(dat,cols=NULL,win=NULL,conv=1,cormethod=1,output=T,pl=1,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Your data frame containing stratigraphic data; any number of columns (variables) are permitted, but the first column should be a location identifier (e.g., depth, height, time).}
  \item{cols}{A vector that identifies the two variable columns to be extracted (first column automatically extracted).}
  \item{win}{Moving window size in units of space or time.}
  \item{conv}{Convention for window placement: (1) center each window on a stratigraphic level in 'dat' (DEFAULT), (2) start with the smallest location datum in 'dat', (3) start with the largest location datum in 'dat'. For options 2 and 3, the center of the window will not necessarily coincide with a measured stratigraphic level in 'dat', but edges of the data set are better preserved.}
  \item{cormethod}{Method used for calculation of correlation coefficient (1=Pearson, 2=Spearman, 3=Kendall)}
  \item{output}{Output results? (T or F)}
  \item{pl}{(1) Plot results at center of window, or (2) create "string of points plot" as in Sageman and Hollander (1999)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\references{
B.B. Sageman and D.H. Hollander, 1999,
\emph{Cross correlation of paleoecological and geochemical proxies: A holistic approach to the study of past global change}, in E. Barrera and C.C. Johnson, eds., GSA Special Paper 332, p. 365-384.
}

\examples{
# generate example series
ex <- cycles(freqs=c(1/40,1/20),noisevar=.2)

# add second variable
ex[3] <- cycles(freqs=c(1/40,1/20),noisevar=0.2)[2]

# jitter sampling times
ex[1]=ex[1]+rnorm(500,sd=1)
# sort
ex = ex[order(ex[,1],na.last=NA,decreasing=FALSE),]

# run mwCor
mwCor(ex,win=50)
}