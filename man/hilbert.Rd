\name{hilbert}
\alias{hilbert}
\title{Hilbert transform of stratigraphic series}
\description{
Calculate instantaneous amplitude (envelope) and phase via Hilbert Transform of stratigraphic series
}
\usage{
hilbert(dat,phase=F,padfac=2,demean=T,detrend=F,output=T,addmean=F,genplot=T,check=T,
        verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series to Hilbert Transform. First column should be location (e.g., depth), second column should be data value.}
  \item{phase}{Calculate instantaneous phase? (T or F)}
  \item{padfac}{Pad with zeros to (padfac*npts) points, where npts is the original number of data points.}
  \item{demean}{Remove mean from data series? (T or F)}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{output}{Return results as new data frame? (T or F)}
  \item{addmean}{Add mean value to instantaneous amplitude? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{check}{Conduct compliance checks before processing? (T or F) In general this should be activated; the option is included for Monte Carlo simulation.}
  \item{verbose}{Verbose output? (T or F)}
}

\examples{
# generate example series with 3 precession terms and noise
ex <- cycles(noisevar=.0004,dt=5)
# bandpass precession terms using cosine-tapered window 
res_ex <- bandpass(ex,flow=0.038,fhigh=0.057,win=2,p=.4)
# hilbert transform
hil_ex <- hilbert(res_ex)
}