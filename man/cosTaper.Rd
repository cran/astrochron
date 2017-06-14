\name{cosTaper}
\alias{cosTaper}
\title{Apply cosine taper to stratigraphic series}
\description{
Apply a "percent-tapered" cosine taper (a.k.a. Tukey window) to a stratigraphic series.
}
\usage{
cosTaper(dat,p=.25,rms=T,demean=T,detrend=F,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for tapering. First column should be location (e.g., depth), second column should be data value. If no data is identified, will output a 256 point taper to evaluate the spectral properties of the window.}
  \item{p}{Cosine-tapered window parameter: p is the percent of the data series tapered (choose 0-1). When p=1, this is equivalent to a Hann taper.}
  \item{rms}{Normalize taper to RMS=1 to preserve power for white process? (T or F)}
  \item{demean}{Remove mean from data series? (T or F)}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{dpssTaper}}, \code{\link{gausTaper}}, and \code{\link{hannTaper}} }