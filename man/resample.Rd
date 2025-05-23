\name{resample}
\alias{resample}
\title{Resample stratigraphic series}
\description{
Resample a stratigraphic series using a new (variably sampled) time or space axis. Values are piecewise-linearly interpolated from original data.
}
\usage{
resample(dat,xout,genplot=T,check=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for resampling. First column should be location (e.g., depth), second column should be data value.}
  \item{xout}{Vector of new sampling locations.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{check}{Conduct compliance checks before processing? (T or F). In general this should be activated; the option is included for Monte Carlo simulation.}
  \item{verbose}{Verbose output? (T or F)}
}