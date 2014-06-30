\name{resample}
\alias{resample}
\title{Resample Stratigraphic Series}
\description{
Resample a stratigraphic series using a new (variably sampled) time or space axis. Values are piecewise-linearly interpolated from original data.
}
\usage{
resample(dat,xout,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for resampling. First column should be location (e.g., depth), second column should be data value.}
  \item{xout}{Vector of new sampling locations.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}