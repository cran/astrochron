\name{detrend}
\alias{detrend}
\title{Subtract linear trend from stratigraphic series}
\description{
Remove linear trend from stratigraphic series
}
\usage{
detrend(dat,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for linear detrending. First column should be location (e.g., depth), second column should be data value.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{arcsinT}}, \code{\link{demean}}, \code{\link{divTrend}}, \code{\link{logT}}, \code{\link{prewhiteAR}}, and \code{\link{prewhiteAR1}}}