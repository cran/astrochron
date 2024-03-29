\name{peak}
\alias{peak}
\title{Identify maxima of peaks in series, filter at desired threshold value}
\description{
Identify maxima of peaks in any 1D or 2D series, filter at desired threshold value.
}
\usage{
peak(dat,level,plateau=F,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{1 or 2 dimensional series. If 2 dimesions, first column should be location (e.g., depth), second column should be data value.}
  \item{level}{Threshold level for filtering peaks. By default all peak maxima reported.}
  \item{plateau}{Output plateau points not evaluated? If T, identified peaks will not be output. (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\examples{
ex=cycles(genplot=FALSE)
peak(ex,level=0.02)
}