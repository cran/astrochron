\name{lowpass}
\alias{lowpass}
\title{Lowpass Filter Stratigraphic Series}
\description{
Lowpass filter stratigraphic series using rectangular, Gaussian or tapered cosine window [cosine window is experimental]
}
\usage{
lowpass(dat,padfac=2,fcut=NULL,win=0,demean=T,detrend=F,alpha=3,p=0.25,xmin=0,xmax=Nyq,
        addmean=T,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for lowpass filtering. First column should be location (e.g., depth), second column should be data value.}
  \item{padfac}{Pad with zeros to (padfac*npts) points, where npts is the original number of data points.}
  \item{fcut}{Cutoff frequency for lowpass filtering.}
  \item{win}{Window type for bandpass filter: 0 = rectangular , 1= Gaussian, 2= Cosine-tapered window.}
  \item{demean}{Remove mean from data series? (T or F)}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{alpha}{Gaussian window parameter: alpha is 1/stdev, a measure of the width of the Dirichlet kernal. Larger values decrease the width of data window, reduce discontinuities, and increase width of the transform. Choose alpha >= 2.5.}
  \item{p}{Cosine-tapered window parameter: p is the percent of the data series tapered (choose 0-1).}
  \item{xmin}{Smallest frequency for plotting.}
  \item{xmax}{Largest frequency for plotting.}
  \item{addmean}{Add mean value to bandpass result? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{bandpass}}, \code{\link{noKernel}}, \code{\link{noLow}}, \code{\link{prewhiteAR}}, and \code{\link{prewhiteAR1}} }
