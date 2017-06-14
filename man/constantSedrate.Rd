\name{constantSedrate}
\alias{constantSedrate}
\title{Apply a constant sedimentation rate model to transform a spatial series to temporal series}
\description{
Apply a constant sedimentation rate model to transform a spatial series to temporal series.
}
\usage{
constantSedrate(dat,sedrate,begin=0,timeDir=1,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series. First column should be location (e.g., depth), second column should be data value.}
  \item{sedrate}{Sedimentation rate, in same spatial units as dat.}
  \item{begin}{Time value to assign to first datum.}
  \item{timeDir}{Direction of floating time in tuned record: 1 = elapsed time increases with depth/height; -1 = elapsed time decreases with depth/height)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}
