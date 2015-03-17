\name{trough}
\alias{trough}
\title{Identify minima of troughs in series, filter at desired threshold value}
\description{
Identify minima of troughs in any 1D or 2D series, filter at desired threshold value.
}
\usage{
trough(dat,level,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{1 or 2 dimensional series. If 2 dimesions, first column should be location (e.g., depth), second column should be data value.}
  \item{level}{Threshold level for filtering troughs. By default all trough minima reported.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\examples{
\donttest{
ex=cycles(genplot=FALSE)
trough(ex,level=-0.02)
}
}