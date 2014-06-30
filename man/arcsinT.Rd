\name{arcsinT}
\alias{arcsinT}
\title{Arcsine Transformation of Stratigraphic Series}
\description{
Arcsine transformation of stratigraphic series
}
\usage{
arcsinT(dat,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for arcsine transformation. First column should be location (e.g., depth), second column should be data value for transformation.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{demean}}, \code{\link{detrend}}, \code{\link{divTrend}}, \code{\link{logT}}, \code{\link{prewhiteAR}}, and \code{\link{prewhiteAR1}}}