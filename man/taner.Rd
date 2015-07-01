\name{taner}
\alias{taner}
\title{Apply Taner Bandpass or Lowpass Filter to Stratigraphic Series}
\description{
Apply Taner bandpass or lowpass filter to stratigraphic series.
}
\usage{
taner(dat,padfac=2,flow=NULL,fhigh=NULL,roll=10^3,demean=T,detrend=F,addmean=T,
       output=1,xmin=0,xmax=Nyq,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for bandpass filtering. First column should be location (e.g., depth), second column should be data value.}
  \item{padfac}{Pad with zeros to (padfac*npts) points, where npts is the original number of data points.}
  \item{flow}{Lowest frequency cut-off (half power point). If this value is not set (NULL), it will default to -1*fhigh, which will create a lowpass filter.}
  \item{fhigh}{Highest frequency cut-off (half power point).}
  \item{roll}{Roll-off rate, in dB/octave.  Typical values are 10^3 to 10^12, but can be larger.}
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

\references{
http://www.rocksolidimages.com/pdf/attrib_revisited.htm#_Toc328470897
}

\seealso{ \code{\link{bandpass}}, \code{\link{lowpass}}, \code{\link{noKernel}}, \code{\link{noLow}}, \code{\link{prewhiteAR}}, and \code{\link{prewhiteAR1}} }

\examples{
\donttest{
# generate example series with 3 precession terms and noise
ex <- cycles(noisevar=.0004,dt=5)
# bandpass precession terms using Taner window 
res_ex <- taner(ex,flow=0.038,fhigh=0.057)

# generate example series with periods of 405 ka, 100 ka and 20 ka, plus noise
ex2=cycles(freqs=c(1/405,1/100,1/20),noisevar=.1,dt=5)
# lowpass filter using Taner window
res_ex2=taner(ex2,fhigh=.02,roll=10^4)
}
}