\name{bandpass}
\alias{bandpass}
\title{Bandpass filter stratigraphic series}
\description{
Bandpass filter stratigraphic series using rectangular, Gaussian or tapered cosine (a.k.a. Tukey) window. This function can also be used to notch fiter a record (see examples).
}
\usage{
bandpass(dat,padfac=2,flow=NULL,fhigh=NULL,win=0,alpha=3,p=0.25,demean=T,
         detrend=F,addmean=T,output=1,xmin=0,xmax=Nyq,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for bandpass filtering. First column should be location (e.g., depth), second column should be data value.}
  \item{padfac}{Pad with zeros to (padfac*npts) points, where npts is the original number of data points.}
  \item{flow}{Lowest frequency to bandpass.}
  \item{fhigh}{Highest frequency to bandpass.}
  \item{win}{Window type for bandpass filter: 0 = rectangular , 1= Gaussian, 2= Cosine-tapered window (a.k.a. Tukey window).}
  \item{alpha}{Gaussian window parameter: alpha is 1/stdev, a measure of the width of the Dirichlet kernel. Choose alpha >= 2.5.}
  \item{p}{Cosine-tapered (Tukey) window parameter: p is the percent of the data series tapered (choose 0-1).}
  \item{demean}{Remove mean from data series? (T or F)}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{addmean}{Add mean value to bandpass result? (T or F)}
  \item{output}{Output: (1) filtered series, (2) bandpass filter window.}
  \item{xmin}{Smallest frequency for plotting.}
  \item{xmax}{Largest frequency for plotting.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}
\value{
bandpassed stratigraphic series.
}

\seealso{ \code{\link{lowpass}}, \code{\link{noKernel}}, \code{\link{noLow}}, \code{\link{prewhiteAR}}, \code{\link{prewhiteAR1}}, and \code{\link{taner}} }

\examples{
# generate example series with periods of 405 ka, 100 ka, 40ka, and 20 ka, plus noise
ex=cycles(freqs=c(1/405,1/100,1/40,1/20),end=1000,dt=5,noisevar=.1)

# bandpass precession term using cosine-tapered window 
bandpass_ex <- bandpass(ex,flow=0.045,fhigh=0.055,win=2,p=.4)

# notch filter (remove) obliquity term using cosine-tapered window
#  if you'd like the final notch filtered record to be centered on the mean proxy 
#  value, set addmean=FALSE
notch_ex <- bandpass(ex,flow=0.02,fhigh=0.03,win=2,p=.4,addmean=FALSE)
notch_ex[2] <- ex[2]-notch_ex[2]
pl(2)
plot(ex,type="l",main="Eccentricity+Obliquity+Precession")
plot(notch_ex,type="l",main="Following application of obliquity notch filter")
}
