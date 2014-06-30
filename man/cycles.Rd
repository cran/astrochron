\name{cycles}
\alias{cycles}
\title{Generate Harmonic Model}
\description{
Make a time series with specified harmonic components and noise
}
\usage{
cycles(freqs=NULL,phase=NULL,amp=NULL,start=0,end=499,dt=1,noisevar=0,genplot=T,verbose=T)
}
\arguments{
  \item{freqs}{Vector with frequencies to model ('linear' frequencies).}
  \item{phase}{Vector with phases for each frequency (phase in radians). Phases are subtracted.}
  \item{amp}{Vector with amplitudes for each frequency.}
  \item{start}{First time/depth/height for output.}
  \item{end}{Last time/depth/height for output.}
  \item{dt}{Sampling interval.}
  \item{noisevar}{Variance of additive Gaussian noise.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\value{
modeled time series.
}

\examples{
## test signal on pg 38 of Choudhury, Shah, and Thornhill (2008)
freqs=c(0.12,0.18,0.30,0.42)
phase=c(-pi/3,-pi/12,-pi/4,-3*pi/8)
amp=c(1,1,1,1)

cycles(freqs,phase,amp,start=0,end=4095,dt=1,noisevar=0.2)
}