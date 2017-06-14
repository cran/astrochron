\name{tune}
\alias{tune}
\title{Tune stratigraphic series}
\description{
Tune stratigraphic series from space to time, using specified control points
}
\usage{
tune(dat,controlPts,extrapolate=F,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for tuning. First column should be location (e.g., depth), second column should be data value.}
  \item{controlPts}{Tuning control points. A data frame or matrix containing two columns: depth, time}
  \item{extrapolate}{Extrapolate sedimentation rates above and below 'tuned' interval? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}