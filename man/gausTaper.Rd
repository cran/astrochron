\name{gausTaper}
\alias{gausTaper}
\title{Apply Gaussian taper to stratigraphic series}
\description{
Apply a Gaussian taper to a stratigraphic series
}
\usage{
gausTaper(dat,alpha=3,rms=T,demean=T,detrend=F,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for tapering. First column should be location (e.g., depth), second column should be data value. If no data is identified, will output a 256 point taper to evaluate the spectral properties of the window.}
  \item{alpha}{Gaussian window parameter: alpha is 1/stdev, a measure of the width of the Dirichlet kernel. Larger values decrease the width of data window, reduce discontinuities, and increase width of the transform. Choose alpha >= 2.5.}
  \item{rms}{Normalize taper to RMS=1 to preserve power for white process? (T or F)}
  \item{demean}{Remove mean from data series? (T or F)}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\references{
Harris, 1978,
\emph{On the use of windows for harmonic analysis with the discrete Fourier transform}: Proceedings of the IEEE, v. 66, p. 51-83.
}

\seealso{ \code{\link{cosTaper}}, \code{\link{dpssTaper}}, and \code{\link{hannTaper}} }