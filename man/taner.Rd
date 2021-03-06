\name{taner}
\alias{taner}
\title{Apply Taner bandpass or lowpass filter to stratigraphic series}
\description{
Apply Taner bandpass or lowpass filter to stratigraphic series. This function can also be used to notch filter or highpass a record (see examples).
}
\usage{
taner(dat,padfac=2,flow=NULL,fhigh=NULL,roll=10^3,demean=T,detrend=F,addmean=T,
       output=1,xmin=0,xmax=Nyq,genplot=T,check=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for bandpass filtering. First column should be location (e.g., depth), second column should be data value.}
  \item{padfac}{Pad with zeros to (padfac*npts) points, where npts is the original number of data points.}
  \item{flow}{Low frequency cut-off for Taner filter (half power point). If this value is not set (NULL), it will default to -1*fhigh, which will create a lowpass filter.}
  \item{fhigh}{High frequency cut-off for Taner filter (half power point).}
  \item{roll}{Roll-off rate, in dB/octave.  Typical values are 10^3 to 10^12, but can be larger.}
  \item{demean}{Remove mean from data series? (T or F)}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{addmean}{Add mean value to bandpass result? (T or F)}
  \item{output}{Output: (1) filtered series, (2) bandpass filter window.}
  \item{xmin}{Smallest frequency for plotting.}
  \item{xmax}{Largest frequency for plotting.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{check}{Conduct compliance checks before processing? (T or F) In general this should be activated; the option is included for Monte Carlo simulation.}
  \item{verbose}{Verbose output? (T or F)}
}

\value{
bandpassed stratigraphic series.
}

\references{
http://www.rocksolidimages.com/pdf/attrib_revisited.htm#_Toc328470897
}

\seealso{ \code{\link{bandpass}}, \code{\link{lowpass}}, \code{\link{noKernel}}, \code{\link{noLow}}, \code{\link{prewhiteAR}}, and \code{\link{prewhiteAR1}} }

\examples{
# generate example series with periods of 405 ka, 100 ka, 40ka, and 20 ka, plus noise
ex=cycles(freqs=c(1/405,1/100,1/40,1/20),end=1000,dt=5,noisevar=.1)

# bandpass precession term using Taner window 
bandpass_ex <- taner(ex,flow=0.045,fhigh=0.055,roll=10^10)

# lowpass filter eccentricity terms using Taner window
lowpass_ex=taner(ex,fhigh=.02,roll=10^10)

# notch filter (remove) obliquity term using Taner window
#  if you'd like the final notch filtered record to be centered on the mean proxy 
#  value, set addmean=FALSE
notch_ex <- taner(ex,flow=0.02,fhigh=0.03,roll=10^10,addmean=FALSE)
notch_ex[2] <- ex[2]-notch_ex[2]
pl(2)
plot(ex,type="l",main="Eccentricity+Obliquity+Precession")
plot(notch_ex,type="l",main="Following application of obliquity notch filter")

# highpass filter obliquity and precession terms using Taner window
#  if you'd like the final highpass filtered record to be centered on the mean proxy 
#  value, set addmean=FALSE
highpass_ex=taner(ex,fhigh=.02,roll=10^10,addmean=FALSE)
highpass_ex[2] <- ex[2]-highpass_ex[2]
pl(2)
plot(ex,type="l",main="Eccentricity+Obliquity+Precession")
plot(highpass_ex,type="l",main="Obliquity+Precession highpassed signal")
}