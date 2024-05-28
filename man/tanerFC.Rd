\name{tanerFC}
\alias{tanerFC}
\title{Apply Taner bandpass or lowpass filter to Fourier coefficients}
\description{
Apply Taner bandpass or lowpass filter to Fourier coefficients.
}
\usage{
tanerFC(fc,npts,flow=NULL,fhigh=NULL,roll=10^3,output=1,genplot=T,verbose=T)
}
\arguments{
  \item{fc}{Fourier coefficients, as output by the function 'periodogram'. The first column is frequency, the second column contains the real coefficients, and the third column contains the imaginary coefficients.}
  \item{npts}{The number of points in the stratigraphic series used to estimate the Fourier coefficients.}
  \item{flow}{Low frequency cut-off for Taner filter (half power point). If this value is not set (NULL), it will default to -1*fhigh, which will create a lowpass filter.}
  \item{fhigh}{High frequency cut-off for Taner filter (half power point).}
  \item{roll}{Roll-off rate, in dB/octave.  Typical values are 10^3 to 10^12, but can be larger.}
  \item{output}{Output: (1) filtered series, (2) bandpass filter window.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
This function is designed for cases when one needs to efficiently apply a range of filter 
parameters to a stratigraphic series. It is used within function 'timeOpt' to speed up 
processing. For more general use, function 'taner' is preferred, which also provides
a range of plots for evaluating the filtering.

Fourier coefficients for 'tanerFC' can be determined with the function 'periodogram',
using options 'output=2' and 'nrm=0'.
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

# calculate the Fourier coefficients using periodogram function. this must be done with 
# no normalization (nrm=0)
fc_ex <- periodogram(ex,demean=TRUE,output=2,nrm=0)

# bandpass precession term using Taner window 
bandpass_ex <- tanerFC(fc=fc_ex,npts=201,flow=0.045,fhigh=0.055,roll=10^10)

# lowpass filter eccentricity terms using Taner window
lowpass_ex <- tanerFC(fc=fc_ex,npts=201,fhigh=.02,roll=10^10)
}