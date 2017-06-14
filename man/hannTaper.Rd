\name{hannTaper}
\alias{hannTaper}
\title{Apply Hann taper to stratigraphic series}
\description{
Apply a Hann (Hanning) taper to a stratigraphic series
}
\usage{
hannTaper(dat,rms=T,demean=T,detrend=F,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for tapering. First column should be location (e.g., depth), second column should be data value. If no data is identified, will output a 256 point taper to evaluate the spectral properties of the window.}
  \item{rms}{Normalize taper to RMS=1 to preserve power for white process? (T or F)}
  \item{demean}{Remove mean from data series? (T or F)}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{cosTaper}}, \code{\link{dpssTaper}}, and \code{\link{gausTaper}} }