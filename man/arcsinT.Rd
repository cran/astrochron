\name{arcsinT}
\alias{arcsinT}
\title{Arcsine transformation of stratigraphic series}
\description{
Arcsine transformation of stratigraphic series
}
\usage{
arcsinT(dat,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for arcsine transformation. Input can have any number of columns desired. If two or more columns are input, the first column must be location (e.g., depth), while remaining columns are data values for transformation.}
  \item{genplot}{Generate summary plots? (T or F). This is automatically deactivated if more than one variable is transformed.}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{demean}}, \code{\link{detrend}}, \code{\link{divTrend}}, \code{\link{logT}}, \code{\link{prewhiteAR}}, and \code{\link{prewhiteAR1}}}