\name{divTrend}
\alias{divTrend}
\title{Divide by Linear Trend in Stratigraphic Series}
\description{
Divide data series value by linear trend observed in stratigraphic series
}
\usage{
divTrend(dat,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for div-trending. First column should be location (e.g., depth), second column should be data value.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{arcsinT}}, \code{\link{demean}}, \code{\link{detrend}}, \code{\link{logT}}, \code{\link{prewhiteAR}}, and \code{\link{prewhiteAR1}}}