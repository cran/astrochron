\name{lowpass}
\alias{lowpass}
\title{Lowpass filter stratigraphic series}
\description{
Lowpass filter stratigraphic series using rectangular, Gaussian or tapered cosine window. This function can also be used to highpass filter a record (see examples).
}
\usage{
lowpass(dat,padfac=2,fcut=NULL,win=0,demean=T,detrend=F,addmean=T,alpha=3,p=0.25,
        xmin=0,xmax=Nyq,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for lowpass filtering. First column should be location (e.g., depth), second column should be data value.}
  \item{padfac}{Pad with zeros to (padfac*npts) points, where npts is the original number of data points.}
  \item{fcut}{Cutoff frequency for lowpass filtering.}
  \item{win}{Window type for bandpass filter: 0 = rectangular , 1= Gaussian, 2= Cosine-tapered window.}
  \item{demean}{Remove mean from data series? (T or F)}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{addmean}{Add mean value to bandpass result? (T or F)}
  \item{alpha}{Gaussian window parameter: alpha is 1/stdev, a measure of the width of the Dirichlet kernal. Larger values decrease the width of data window, reduce discontinuities, and increase width of the transform. Choose alpha >= 2.5.}
  \item{p}{Cosine-tapered window parameter: p is the percent of the data series tapered (choose 0-1).}
  \item{xmin}{Smallest frequency for plotting.}
  \item{xmax}{Largest frequency for plotting.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{bandpass}}, \code{\link{noKernel}}, \code{\link{noLow}}, \code{\link{prewhiteAR}}, \code{\link{prewhiteAR1}}, and \code{\link{taner}} }

\examples{
# generate example series with periods of 405 ka, 100 ka, 40ka, and 20 ka, plus noise
ex=cycles(freqs=c(1/405,1/100,1/40,1/20),end=1000,dt=5,noisevar=.1)

# lowpass filter eccentricity terms using cosine-tapered window
lowpass_ex=lowpass(ex,fcut=.02,win=2,p=.4)

# highpass filter obliquity and precession terms using cosine-tapered window
#  if you'd like the final notch filtered record to be centered on the mean proxy 
#  value, set addmean=FALSE
highpass_ex=lowpass(ex,fcut=.02,win=2,p=.4,addmean=FALSE)
highpass_ex[2] <- ex[2]-highpass_ex[2]
pl(2)
plot(ex,type="l",main="Eccentricity+Obliquity+Precession")
plot(highpass_ex,type="l",main="Obliquity+Precession highpassed signal")

}