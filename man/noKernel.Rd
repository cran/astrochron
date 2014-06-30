\name{noKernel}
\alias{noKernel}
\title{Remove Gaussian Kernel Smoother from Stratigraphic Series}
\description{
Estimate trend and remove from stratigraphic series using a Gaussian kernel smoother }
\usage{
noKernel(dat,smooth=0.1,sort=F,output=1,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for smoothing. First column should be location (e.g., depth), second column should be data value.}
  \item{smooth}{Degree of smoothing with a Gaussian kernal (0 = no smoothing); for a value of 0.5, the kernel is scaled so that its quartiles (viewed as prob densities) are at +/- 25 percent of the data series length. Must be > 0.}
  \item{sort}{Sort data into increasing depth (required for ksmooth)? (T or F)}
  \item{output}{1= output residual values; 2= output Gaussian kernel smoother.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{bandpass}}, \code{\link{lowpass}}, \code{\link{noLow}}, \code{\link{prewhiteAR}}, and \code{\link{prewhiteAR1}} }
