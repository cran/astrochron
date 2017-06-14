\name{flip}
\alias{flip}
\title{Flip stratigraphic series}
\description{
Flip the stratigraphic order of your data series (e.g., convert stratigraphic depth series to 
height series, relative to a defined datum.)
}
\usage{
flip(dat,begin=0,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series. First column should be location (e.g., depth), second column should be data value.}
  \item{begin}{Depth/height value to assign to (new) first stratigraphic datum.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}
