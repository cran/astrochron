\name{s}
\alias{s}
\title{Standardize variable in Stratigraphic Series}
\description{
Standardize variable in Stratigraphic Series (subtract mean value and divide by standard deviation)
}
\usage{
s(dat,genplot=F,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for standardization. First column should be location (e.g., depth), second column should be data value.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}