\name{ar1}
\alias{ar1}
\title{Generate AR(1) model}
\description{
Generate an AR(1) model
}
\usage{
ar1(npts=1024,dt=1,mean=0,sdev=1,rho=0.9,genplot,verbose)

}
\arguments{
  \item{npts}{number of time series data points}
  \item{dt}{sampling interval}
  \item{mean}{mean value for Gaussian deviates}
  \item{sdev}{standard deviation for Gaussian deviates}
  \item{rho}{AR(1) coefficient}
  \item{genplot}{generate summary plots (T or F)}
  \item{verbose}{verbose output (T or F)}
}

\examples{
# generate an AR1 series with 500 points, sampling interval of 5 ka, and rho of 0.8
ar1(npts=500,dt=5,rho=0.8)
}
