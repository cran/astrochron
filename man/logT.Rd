\name{logT}
\alias{logT}
\title{Log Transformation of Stratigraphic Series}
\description{
Log transformation of stratigraphic series.
}
\usage{
logT(dat,c=0,opt=1,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for log transformation. First column should be location (e.g., depth), second column should be data value for transformation.}
  \item{c}{Constant to add prior to log transformation. Default = 0.}
  \item{opt}{(1) use natural logarithm, (2) use log10. Default = 1.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{arcsinT}}, \code{\link{demean}}, \code{\link{detrend}}, \code{\link{divTrend}}, \code{\link{prewhiteAR}}, and \code{\link{prewhiteAR1}}}