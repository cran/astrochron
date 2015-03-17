\name{hilbert}
\alias{hilbert}
\title{Hilbert Transform of Stratigraphic Series}
\description{
Calculate instantaneous amplitude (envelope) via Hilbert Transform of stratigraphic series
}
\usage{
hilbert(dat,padfac=2,demean=T,detrend=F,output=T,addmean=F,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series to Hilbert Transform. First column should be location (e.g., depth), second column should be data value.}
  \item{padfac}{Pad with zeros to (padfac*npts) points, where npts is the original number of data points.}
  \item{demean}{Remove mean from data series? (T or F)}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{output}{Return results as new data frame? (T or F)}
  \item{addmean}{Add mean value to instantaneous amplitude? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\examples{
\donttest{
# generate example series with 3 precession terms and noise
ex <- cycles(noisevar=.0004,dt=5)
# bandpass precession terms using cosine-tapered window 
res_ex <- bandpass(ex,flow=0.038,fhigh=0.057,win=2,p=.4)
# hilbert transform
hil_ex <- hilbert(res_ex)
}
}