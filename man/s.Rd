\name{s}
\alias{s}
\title{Standardize variable in stratigraphic series}
\description{
Standardize variable in stratigraphic series (subtract mean value and divide by standard deviation)
}
\usage{
s(dat,genplot=F,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for standardization. First column should be location (e.g., depth), second column should be data value.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}