\name{prewhiteAR1}
\alias{prewhiteAR1}
\title{Prewhiten Stratigraphic Series with AR1 filter, using 'Standard' or Unbiased Estimate of rho}
\description{
Prewhiten stratigraphic series using autoregressive-1 (AR1) filter. Rho can be estimated using the 'standard' approach, or following a bias correction.
}
\usage{
prewhiteAR1(dat,setrho=NULL,bias=F,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for prewhitening. First column should be location (e.g., depth), second column should be data value for prewhitening. Series must have uniform sampling interval.}
  \item{setrho}{Specified lag-1 correlation coefficient (rho). By default, rho is calculated.}
  \item{bias}{Calculate unbiased estimate of rho, as in Mudelsee (2010, eq. 2.45). (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}
\references{
M. Mudelsee, 2010,
\emph{Climate Time Series Analysis: Classical Statistical and Bootstrap Methods}, 474 pp., Springer, Dordrecht, Netherlands.
}

\seealso{ \code{\link{arcsinT}}, \code{\link{bandpass}}, \code{\link{demean}}, \code{\link{detrend}}, \code{\link{divTrend}}, \code{\link{logT}}, \code{\link{lowpass}}, \code{\link{noKernel}}, and \code{\link{prewhiteAR}} }