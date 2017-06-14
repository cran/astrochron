\name{logT}
\alias{logT}
\title{Log transformation of stratigraphic series}
\description{
Log transformation of stratigraphic series.
}
\usage{
logT(dat,c=0,opt=1,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for log transformation. Input can have any number of columns desired. If two or more columns are input, the first column must be location (e.g., depth), while remaining columns are data values for transformation.}
  \item{c}{Constant to add prior to log transformation. Default = 0.}
  \item{opt}{(1) use natural logarithm, (2) use log10. Default = 1.}
  \item{genplot}{Generate summary plots? (T or F). This is automatically deactivated if more than one variable is transformed.}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{arcsinT}}, \code{\link{demean}}, \code{\link{detrend}}, \code{\link{divTrend}}, \code{\link{prewhiteAR}}, and \code{\link{prewhiteAR1}}}