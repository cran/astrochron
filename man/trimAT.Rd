\name{trimAT}
\alias{trimAT}
\title{Remove Outliers from Stratigraphic Series}
\description{
Remove outliers from stratigraphic series, using specified threshold value.
}
\usage{
trimAT(dat,thresh=0,dir=2,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for outlier removal. First column should be location (e.g., depth), second column should be data value.}
  \item{thresh}{Threshold value for outlier detection.}
  \item{dir}{Remove values (1) smaller than or (2) larger than this threshold?}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{delPts}}, \code{\link{idPts}}, \code{\link{iso}} and \code{\link{trim}} }