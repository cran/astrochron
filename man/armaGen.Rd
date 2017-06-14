\name{armaGen}
\alias{armaGen}
\title{Generate autoregressive moving-average model}
\description{
Generate an autoregressive moving-average time series model
}
\usage{
armaGen(npts=1024,dt=1,m=0,std=1,rhos=c(0.9),thetas=c(0),genplot=T,verbose=T)

}
\arguments{
  \item{npts}{Number of time series data points.}
  \item{dt}{Sampling interval.}
  \item{m}{Mean value of final time series.}
  \item{std}{Standard deviation of final time series.}
  \item{rhos}{Vector of AR coefficients for each order.}
  \item{thetas}{Vector of MA coefficients for each order.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}