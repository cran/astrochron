\name{dpssTaper}
\alias{dpssTaper}
\title{Apply DPSS taper to stratigraphic series}
\description{
Apply a single Discrete Prolate Spheroidal Sequence (DPSS) taper to a stratigraphic series
}
\usage{
dpssTaper(dat,tbw=1,num=1,rms=T,demean=T,detrend=F,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for tapering. First column should be location (e.g., depth), second column should be data value. If no data is identified, will output a 256 point taper to evaluate the spectral properties of the window.}
  \item{tbw}{Time-bandwidth product for the DPSS}
  \item{num}{Which one of the DPSS would you like to use?}
  \item{rms}{Normalize taper to RMS=1 to preserve power for white process? (T or F)}
  \item{demean}{Remove mean from data series? (T or F)}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{cosTaper}}, \code{\link{gausTaper}}, and \code{\link{hannTaper}} }