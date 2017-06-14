\name{noLow}
\alias{noLow}
\title{Fit and remove Lowess smoother from stratigraphic series}
\description{
Fit and remove lowess smoother from stratigraphic series
}
\usage{
noLow(dat,smooth=.20,output=1,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for lowess smoother removal. First column should be location (e.g., depth), second column should be data value.}
  \item{smooth}{Lowess smoothing parameter.}
  \item{output}{1= output residual values; 2= output lowess fit}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{bandpass}}, \code{\link{lowpass}}, \code{\link{noKernel}}, \code{\link{prewhiteAR}}, and \code{\link{prewhiteAR1}} }
