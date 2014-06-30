\name{bandpass}
\alias{bandpass}
\title{Bandpass Filter Stratigraphic Series}
\description{
Bandpass filter stratigraphic series using rectangular, Gaussian or tapered cosine window
}
\usage{
bandpass(dat,padfac=2,flow=NULL,fhigh=NULL,win=0,alpha=3,p=0.25,demean=T,
         detrend=F,xmin=0,xmax=Nyq,addmean=T,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for bandpass filtering. First column should be location (e.g., depth), second column should be data value.}
  \item{padfac}{Pad with zeros to (padfac*npts) points, where npts is the original number of data points.}
  \item{flow}{Lowest frequency to bandpass.}
  \item{fhigh}{Highest frequency to bandpass.}
  \item{win}{Window type for bandpass filter: 0 = rectangular , 1= Gaussian, 2= Cosine-tapered window.}
  \item{alpha}{Gaussian window parameter: alpha is 1/stdev, a measure of the width of the Dirichlet kernel. Choose alpha >= 2.5.}
  \item{p}{Cosine-tapered window parameter: p is the percent of the data series tapered (choose 0-1).}
  \item{demean}{Remove mean from data series? (T or F)}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{xmin}{Smallest frequency for plotting.}
  \item{xmax}{Largest frequency for plotting.}
  \item{addmean}{Add mean value to bandpass result? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}
\value{
bandpassed stratigraphic series.
}

\seealso{ \code{\link{lowpass}}, \code{\link{noKernel}}, \code{\link{noLow}}, \code{\link{prewhiteAR}}, and \code{\link{prewhiteAR1}} }

\examples{
# generate example series with 3 precession terms and noise
ex <- cycles(noisevar=.0004,dt=5)
# bandpass precession terms using cosine-tapered window 
res_ex <- bandpass(ex,flow=0.038,fhigh=0.057,win=2,p=.4)
}
