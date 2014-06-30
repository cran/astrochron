\name{demean}
\alias{demean}
\title{Remove Mean Value from Stratigraphic Series}
\description{
Remove mean value from stratigraphic series
}
\usage{
demean(dat,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for mean removal. First column should be location (e.g., depth), second column should be data value.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{arcsinT}}, \code{\link{detrend}}, \code{\link{divTrend}}, \code{\link{logT}}, \code{\link{prewhiteAR}}, and \code{\link{prewhiteAR1}}}