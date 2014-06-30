\name{logT}
\alias{logT}
\title{Log Transformation of Stratigraphic Series}
\description{
Log transformation of stratigraphic series.
}
\usage{
logT(dat,c=0,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for log transformation. First column should be location (e.g., depth), second column should be data value for transformation.}
  \item{c}{Constant to add prior to log transformation. Default = 0.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{arcsinT}}, \code{\link{demean}}, \code{\link{detrend}}, \code{\link{divTrend}}, \code{\link{prewhiteAR}}, and \code{\link{prewhiteAR1}}}