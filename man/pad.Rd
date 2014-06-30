\name{pad}
\alias{pad}
\title{Pad Stratigraphic Series with Zeros}
\description{
Pad Stratigraphic Series with Zeros ("zero padding")
}
\usage{
pad(dat,zeros,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for mean removal. First column should be location (e.g., depth), second column should be data value.}
  \item{zeros}{Number of zeros to add on the end of the series. By default, the number of points will be doubled.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}