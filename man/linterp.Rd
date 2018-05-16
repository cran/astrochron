\name{linterp}
\alias{linterp}
\title{Piecewise linear interpolation of stratigraphic series}
\description{
Interpolate stratigraphic series onto a evenly sampled grid, using piecewise linear interpolation
}
\usage{
linterp(dat,dt,start,genplot=T,check=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for piecewise linear interpolation. First column should be location (e.g., depth), second column should be data value.}
  \item{dt}{New sampling interval.}
  \item{start}{Start interpolating at what time/depth/height value? By default, the first value of the stratigraphic series will be used.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{check}{Conduct compliance checks before processing? (T or F) In general this should be activated; the option is included for Monte Carlo simulation.}
  \item{verbose}{Verbose output? (T or F)}
}